"""
RELIO ENGINE
============
Refactored LitMap V0.1 core
- FAST mode: abstracts only
- FULL mode: PMC full text + maze graph + contradictions
- Outputs:
  - python-runner/results/litmap_report_<job>.txt
  - python-runner/results/images/litmap_maze_<job>.png
"""

import os
import time
import re
import math
import requests
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor

from Bio import Entrez

# =========================
# CONFIG
# =========================
Entrez.email = "relio@paperoservices.ai"
Entrez.tool = "Relio"
NCBI_DELAY = 0.35


# =========================
# UTILS
# =========================
def sleep():
    time.sleep(NCBI_DELAY)


def get_context_direction(text, gene):
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts:
        return "neutral"

    idx = starts[0]
    snippet = text[max(0, idx - 120): min(len(text), idx + 120)].lower()

    if any(w in snippet for w in ["increase", "induce", "upregulat", "activat", "stimulat", "enhance", "promot"]):
        return "up"
    if any(w in snippet for w in ["decrease", "inhibit", "suppress", "reduc", "block", "attenuate", "prevent"]):
        return "down"
    return "neutral"


def analyze_contradictions(gene_evidence):
    conflicts = []
    for gene, hits in gene_evidence.items():
        dirs = [h["direction"] for h in hits]
        up = dirs.count("up")
        down = dirs.count("down")
        total = up + down
        if total > 1 and min(up, down) / total > 0.2:
            conflicts.append({
                "gene": gene,
                "up": up,
                "down": down,
                "consensus": "Contradictory"
            })
    return conflicts


def compute_fingerprint(gene_evidence, pathway_map):
    up = sum(1 for hits in gene_evidence.values() for h in hits if h["direction"] == "up")
    down = sum(1 for hits in gene_evidence.values() for h in hits if h["direction"] == "down")

    if up > down * 1.5:
        direction = "Predominant Activation"
    elif down > up * 1.5:
        direction = "Predominant Inhibition"
    else:
        direction = "Balanced Modulation"

    path_sizes = [len(v) for v in pathway_map.values()]
    total = sum(path_sizes)
    entropy = 0.0
    if total:
        probs = [p / total for p in path_sizes]
        entropy = -sum(p * math.log2(p) for p in probs)

    return {"direction": direction, "entropy": round(entropy, 3)}


# =========================
# GENE WHITELIST
# =========================
def build_gene_whitelist(outcome):
    genes = set()
    try:
        sleep()
        h = Entrez.esearch(
            db="gene",
            term=f"{outcome} AND Homo sapiens[Organism]",
            retmax=200
        )
        ids = Entrez.read(h)["IdList"]
        h.close()

        if not ids:
            return genes

        sleep()
        h = Entrez.esummary(db="gene", id=",".join(ids))
        data = Entrez.read(h, validate=False)
        h.close()

        for d in data["DocumentSummarySet"]["DocumentSummary"]:
            g = d.get("Name", "")
            if g.isupper() and 2 <= len(g) <= 10:
                genes.add(g)

    except Exception:
        pass

    return genes


# =========================
# MINING
# =========================
def mine_abstracts_fast(compound, outcome, genes, limit=100):
    gene_evidence = defaultdict(list)
    pmids = []

    try:
        h = Entrez.esearch(
            db="pubmed",
            term=f"{compound} AND {outcome}",
            retmax=limit,
            sort="relevance"
        )
        pmids = Entrez.read(h)["IdList"]
        h.close()

        if not pmids:
            return {}, []

        sleep()
        h = Entrez.efetch(
            db="pubmed",
            id=",".join(pmids),
            rettype="abstract",
            retmode="text"
        )
        blob = h.read()
        h.close()

        abstracts = [a for a in blob.split("\n\n") if len(a) > 50]

        for i, txt in enumerate(abstracts):
            pid = pmids[i] if i < len(pmids) else "Unknown"
            for g in genes:
                if re.search(rf"\b{g}\b", txt):
                    gene_evidence[g].append({
                        "id": pid,
                        "direction": get_context_direction(txt, g)
                    })

    except Exception:
        pass

    return gene_evidence, pmids


def mine_pmc_full(compound, outcome, genes, limit=100):
    gene_evidence = defaultdict(list)
    pmc_ids = []

    try:
        h = Entrez.esearch(
            db="pmc",
            term=f"{compound} AND {outcome} AND open access[filter]",
            retmax=limit,
            sort="relevance"
        )
        pmc_ids = Entrez.read(h)["IdList"]
        h.close()

        if not pmc_ids:
            return {}, []

        batch = 20
        for i in range(0, len(pmc_ids), batch):
            sleep()
            h = Entrez.efetch(
                db="pmc",
                id=",".join(pmc_ids[i:i + batch]),
                retmode="xml"
            )
            articles = Entrez.read(h, validate=False)
            h.close()

            art_list = articles if isinstance(articles, list) else articles.get("pmc-articles", [])
            if isinstance(art_list, dict):
                art_list = [art_list]

            for art in art_list:
                text = str(art)
                pmc = "Unknown"
                try:
                    for aid in art.get("front", {}).get("article-meta", {}).get("article-id", []):
                        if aid.attributes.get("pub-id-type") == "pmc":
                            pmc = "PMC" + str(aid)
                except Exception:
                    pass

                for g in genes:
                    if re.search(rf"[^a-zA-Z]{g}[^a-zA-Z]", text):
                        gene_evidence[g].append({
                            "id": pmc,
                            "direction": get_context_direction(text, g)
                        })

    except Exception:
        pass

    return gene_evidence, pmc_ids


# =========================
# PATHWAYS
# =========================
def fetch_live_pathways(gene):
    paths = set()
    try:
        r = requests.get(
            "https://reactome.org/ContentService/search/query",
            params={"query": gene, "species": "Homo sapiens", "types": "Pathway"},
            timeout=2
        )
        if r.ok:
            for res in r.json().get("results", []):
                paths.add(res["name"])
    except Exception:
        pass

    try:
        r = requests.get(
            "https://webservice.wikipathways.org/findPathwaysByText",
            params={"query": gene, "format": "json"},
            timeout=2
        )
        if r.ok:
            for res in r.json().get("result", []):
                if res.get("species") == "Homo sapiens":
                    paths.add(res["name"])
    except Exception:
        pass

    return list(paths)


def build_pathway_map(gene_evidence):
    pathway_map = defaultdict(set)
    genes = list(gene_evidence.keys())

    with ThreadPoolExecutor(max_workers=10) as ex:
        results = ex.map(fetch_live_pathways, genes)
        for g, paths in zip(genes, results):
            for p in paths:
                p_clean = p.split(" - ")[0].strip()
                if "disease" not in p_clean.lower():
                    pathway_map[p_clean].add(g)

    return pathway_map


# =========================
# GRAPH
# =========================
def draw_graph(compound, gene_evidence, pathway_map, mode_name, job_id):
    if "FAST" in mode_name:
        return None

    G = nx.Graph()

    if pathway_map:
        top = sorted(pathway_map.items(), key=lambda x: len(x[1]), reverse=True)[:20]
        for p, genes in top:
            for g in genes:
                G.add_edge(g, p)

    else:
        top_genes = sorted(gene_evidence.keys(), key=lambda g: len(gene_evidence[g]), reverse=True)[:20]
        for g in top_genes:
            G.add_edge(compound, g)

    plt.figure(figsize=(18, 14))
    pos = nx.kamada_kawai_layout(G)

    nx.draw(G, pos, with_labels=True, node_size=1200, font_size=8)
    plt.title(f"Relio LitMap: {compound}")

    out = f"python-runner/results/images/litmap_maze_{job_id}.png"
    plt.savefig(out, dpi=300)
    plt.close()

    return out


# =========================
# REPORT
# =========================
def generate_report(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name, job_id):
    fingerprint = compute_fingerprint(gene_evidence, pathway_map)
    conflicts = analyze_contradictions(gene_evidence) if "FULL" in mode_name else []

    lines = []
    lines.append("=" * 80)
    lines.append(f"RELIO LITMAP REPORT: {compound.upper()} + {outcome.upper()}")
    lines.append("=" * 80)
    lines.append(f"Mode          : {mode_name}")
    lines.append(f"Papers        : {len(paper_ids)}")
    lines.append(f"Direction     : {fingerprint['direction']}")
    lines.append("")

    for g, ev in sorted(gene_evidence.items(), key=lambda x: len(x[1]), reverse=True):
        dirs = Counter([e["direction"] for e in ev]).most_common(1)[0][0]
        lines.append(f"{g:8} | {len(ev):3} | {dirs}")

    out = f"python-runner/results/litmap_report_{job_id}.txt"
    with open(out, "w") as f:
        f.write("\n".join(lines))

    return out


# =========================
# MAIN ENTRY
# =========================
def run_litmap(compound, outcome, mode, job_id):
    os.makedirs("python-runner/results/images", exist_ok=True)

    whitelist = build_gene_whitelist(outcome)
    if not whitelist:
        return {"error": "Gene whitelist failed"}

    if mode == "FAST":
        gene_ev, ids = mine_abstracts_fast(compound, outcome, whitelist)
        mode_name = "FAST (Abstracts)"
    else:
        gene_ev, ids = mine_pmc_full(compound, outcome, whitelist)
        mode_name = "FULL (PMC Full Text)"

    if not gene_ev:
        return {"error": "No evidence found"}

    pmap = build_pathway_map(gene_ev)
    graph_path = draw_graph(compound, gene_ev, pmap, mode_name, job_id)
    report_path = generate_report(compound, outcome, gene_ev, pmap, ids, mode_name, job_id)

    return {
        "compound": compound,
        "outcome": outcome,
        "mode": mode_name,
        "genes": len(gene_ev),
        "papers": len(ids),
        "graph": graph_path,
        "report": report_path
    }