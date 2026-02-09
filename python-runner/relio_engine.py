import os
import time
import re
import math
import json
import requests
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez

# =========================
# CONFIG
# =========================
Entrez.email = "your.email@example.com"
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
    snippet = text[max(0, idx-120): min(len(text), idx+120)].lower()

    if any(w in snippet for w in ["increase", "induce", "upregulat", "activat", "stimulat", "enhance", "promot"]):
        return "up"
    if any(w in snippet for w in ["decrease", "inhibit", "suppress", "reduc", "block", "attenuate", "prevent"]):
        return "down"
    return "neutral"


# =========================
# SCIENCE LAYER
# =========================
def compute_fingerprint(gene_evidence):
    ups = sum(1 for hits in gene_evidence.values() for h in hits if h["direction"] == "up")
    downs = sum(1 for hits in gene_evidence.values() for h in hits if h["direction"] == "down")

    total = ups + downs
    if total == 0:
        return {"direction": "Neutral", "entropy": 0, "up": 0, "down": 0}

    p_up = ups / total
    p_down = downs / total
    entropy = -sum(p * math.log2(p) for p in [p_up, p_down] if p > 0)

    if p_up > 0.6:
        direction = "Predominant Activation"
    elif p_down > 0.6:
        direction = "Predominant Inhibition"
    else:
        direction = "Mixed / Uncertain"

    return {
        "direction": direction,
        "entropy": round(entropy, 3),
        "up": ups,
        "down": downs
    }


def analyze_contradictions(gene_evidence):
    contradictions = []
    for gene, hits in gene_evidence.items():
        ups = sum(1 for h in hits if h["direction"] == "up")
        downs = sum(1 for h in hits if h["direction"] == "down")
        if ups > 0 and downs > 0:
            contradictions.append({
                "gene": gene,
                "up": ups,
                "down": downs,
                "status": "Contradictory"
            })
    return contradictions


# =========================
# 1. GENE WHITELIST
# =========================
def build_gene_whitelist(outcome):
    print(f"[1] Building whitelist for {outcome}")
    genes = set()
    try:
        sleep()
        h = Entrez.esearch(db="gene", term=f"{outcome} AND Homo sapiens[Organism]", retmax=200)
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

    except Exception as e:
        print("Whitelist failed:", e)

    print("Whitelist size:", len(genes))
    return genes


# =========================
# 2. ABSTRACT MINING
# =========================
def mine_abstracts_fast(compound, outcome, genes, limit=100):
    gene_evidence = defaultdict(list)
    pmids = []

    try:
        h = Entrez.esearch(db="pubmed", term=f"{compound} AND {outcome}", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        pmids = ids

        if not ids:
            return gene_evidence, pmids

        sleep()
        h = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
        text_blob = h.read()
        h.close()

        abstracts = [a for a in text_blob.split("\n\n") if len(a) > 50]

        for i, txt in enumerate(abstracts):
            pid = ids[i] if i < len(ids) else "Unknown"
            for g in genes:
                if re.search(rf"\b{g}\b", txt):
                    direction = get_context_direction(txt, g)
                    gene_evidence[g].append({"pmid": pid, "direction": direction})

    except Exception as e:
        print("Mining error:", e)

    return gene_evidence, pmids


# =========================
# 3. PATHWAYS
# =========================
def fetch_live_pathways(gene):
    paths = set()
    try:
        r = requests.get("https://reactome.org/ContentService/search/query",
                         params={"query": gene, "species": "Homo sapiens", "types": "Pathway"}, timeout=2)
        if r.ok:
            for res in r.json().get("results", []):
                paths.add(res["name"])
    except:
        pass
    return list(paths)


def build_pathway_map(gene_evidence):
    pathway_map = defaultdict(set)
    genes = list(gene_evidence.keys())

    with ThreadPoolExecutor(max_workers=10) as executor:
        results = executor.map(fetch_live_pathways, genes)
        for gene, paths in zip(genes, results):
            for p in paths:
                p = p.split(" - ")[0]
                if "disease" not in p.lower():
                    pathway_map[p].add(gene)

    return pathway_map


# =========================
# 4. GRAPH
# =========================
def draw_graph(gene_evidence, pathway_map, graph_path):
    G = nx.Graph()

    for p, genes in list(pathway_map.items())[:15]:
        for g in genes:
            G.add_edge(g, p)

    if not G.nodes:
        for g in list(gene_evidence.keys())[:15]:
            G.add_node(g)

    plt.figure(figsize=(14, 12))
    pos = nx.kamada_kawai_layout(G)
    nx.draw(G, pos, with_labels=True, node_size=1400, font_size=9)
    plt.savefig(graph_path, dpi=300)
    plt.close()


# =========================
# MAIN ENGINE
# =========================
def run_relio_job(compound, outcome, mode, job_id):

    whitelist = build_gene_whitelist(outcome)
    gene_evidence, pmids = mine_abstracts_fast(compound, outcome, whitelist)
    pathway_map = build_pathway_map(gene_evidence)

    # Gene stats
    gene_stats = {}
    for g, hits in gene_evidence.items():
        ups = sum(1 for h in hits if h["direction"] == "up")
        downs = sum(1 for h in hits if h["direction"] == "down")
        neutral = sum(1 for h in hits if h["direction"] == "neutral")

        if ups > downs:
            dom = "UP"
        elif downs > ups:
            dom = "DOWN"
        elif ups == downs and ups > 0:
            dom = "MIXED"
        else:
            dom = "NEUTRAL"

        gene_stats[g] = {
            "up": ups,
            "down": downs,
            "neutral": neutral,
            "dominant": dom,
            "pathways": [p for p, gs in pathway_map.items() if g in gs]
        }

    fingerprint = compute_fingerprint(gene_evidence)
    contradictions = analyze_contradictions(gene_evidence)

    os.makedirs("python-runner/results", exist_ok=True)
    os.makedirs("python-runner/images", exist_ok=True)

    graph_file = f"{job_id}.png"
    draw_graph(gene_evidence, pathway_map, f"python-runner/images/{graph_file}")

    result = {
        "job_id": job_id,
        "compound": compound,
        "outcome": outcome,
        "fingerprint": fingerprint,
        "contradictions": contradictions,
        "genes": gene_stats,
        "pathways": {k: list(v) for k, v in pathway_map.items()},
        "sources": [f"https://pubmed.ncbi.nlm.nih.gov/{p}/" for p in pmids],
        "graph": graph_file
    }

    with open(f"python-runner/results/{job_id}.json", "w") as f:
        json.dump(result, f, indent=2)

    return result
