import os
import re
import math
import time
import json
import requests
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict, Counter
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez

# =========================
# Configuration
# =========================
Entrez.email = "your.email@example.com"
Entrez.tool = "Relio_V1"
NCBI_DELAY = 0.35

# =========================
# Utilities
# =========================
def sleep():
    time.sleep(NCBI_DELAY)

def get_context_direction(text, gene):
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts:
        return "neutral"
    idx = starts[0]
    snippet = text[max(0, idx-100): idx+100].lower()
    if any(w in snippet for w in ["increase", "induce", "upregulat", "activat", "stimulat", "enhance", "promot"]):
        return "up"
    if any(w in snippet for w in ["decrease", "inhibit", "suppress", "reduc", "block", "attenuate", "prevent"]):
        return "down"
    return "neutral"

def compute_fingerprint(gene_evidence, pathway_map):
    up = sum(1 for hits in gene_evidence.values() for h in hits if h['direction']=='up')
    down = sum(1 for hits in gene_evidence.values() for h in hits if h['direction']=='down')

    if up > down * 1.5:
        direction = "Predominant Activation"
    elif down > up * 1.5:
        direction = "Predominant Inhibition"
    else:
        direction = "Balanced Modulation"

    path_counts = [len(genes) for genes in pathway_map.values()]
    total_conn = sum(path_counts)
    entropy = 0.0
    if total_conn > 0:
        probs = [c / total_conn for c in path_counts]
        entropy = -sum(p * math.log2(p) for p in probs)

    return {"direction": direction, "entropy": entropy}

# =========================
# Whitelist
# =========================
def build_gene_whitelist(outcome):
    print(f"[1] 🧠 Building Gene Whitelist for '{outcome}'...")
    try:
        sleep()
        h = Entrez.esearch(db="gene", term=f"{outcome} AND Homo sapiens[Organism]", retmax=200)
        ids = Entrez.read(h)["IdList"]
        h.close()
        genes = set()
        if ids:
            sleep()
            h = Entrez.esummary(db="gene", id=",".join(ids))
            data = Entrez.read(h, validate=False)
            h.close()
            for d in data["DocumentSummarySet"]["DocumentSummary"]:
                g = d.get("Name", "")
                if g.isupper() and 2 <= len(g) <= 10:
                    genes.add(g)
        print(f"    → Whitelist: {len(genes)} verified genes.")
        return genes
    except Exception:
        return set()

# =========================
# Mining Engines
# =========================
def mine_abstracts_fast(compound, outcome, genes, limit=100):
    print(f"\n[2] ⚡ FAST MODE: Fetching up to {limit} PubMed Abstracts...")
    gene_evidence = defaultdict(list)
    pmids = []
    try:
        h = Entrez.esearch(db="pubmed", term=f"{compound} AND {outcome}", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        pmids = ids
        if not ids:
            return {}, []
        print(f"    → Found {len(ids)} abstracts. Downloading...")

        sleep()
        h = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
        text_blob = h.read()
        h.close()
        abstracts = [a.strip() for a in text_blob.split("\n\n") if len(a) > 50]

        for i, txt in enumerate(abstracts):
            pid = ids[i] if i < len(ids) else "Unknown"
            for g in genes:
                if re.search(rf"\b{g}\b", txt):
                    direction = get_context_direction(txt, g)
                    gene_evidence[g].append({"id": pid, "direction": direction})
        return gene_evidence, pmids
    except Exception:
        return {}, []

def mine_pmc_full(compound, outcome, genes, limit=100):
    print(f"\n[2] 🤿 FULL MODE: Searching up to {limit} PMC Full Text Articles...")
    gene_evidence = defaultdict(list)
    pmc_list = []
    try:
        h = Entrez.esearch(db="pmc", term=f"{compound} AND {outcome} AND open access[filter]", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        pmc_list = ids
        if not ids:
            return {}, []

        batch_size = 20
        for i in range(0, len(ids), batch_size):
            try:
                sleep()
                h = Entrez.efetch(db="pmc", id=",".join(ids[i:i+batch_size]), retmode="xml")
                articles = Entrez.read(h, validate=False)
                h.close()
                art_list = articles if isinstance(articles, list) else articles.get("pmc-articles", [])
                if isinstance(art_list, dict): art_list = [art_list]
                for art in art_list:
                    current_pmc = "Unknown"
                    try:
                        for aid in art.get('front', {}).get('article-meta', {}).get('article-id', []):
                            if aid.attributes.get('pub-id-type') == 'pmc':
                                current_pmc = "PMC" + str(aid)
                    except: pass
                    full_text = str(art)
                    for g in genes:
                        if re.search(rf"[^a-zA-Z]{g}[^a-zA-Z]", full_text):
                            direction = get_context_direction(full_text, g)
                            gene_evidence[g].append({"id": current_pmc, "direction": direction})
            except: continue
        return gene_evidence, pmc_list
    except: return {}, []

# =========================
# Pathway Mapping
# =========================
def fetch_live_pathways(gene):
    paths = set()
    try:
        r = requests.get("https://reactome.org/ContentService/search/query",
                         params={"query": gene, "species": "Homo sapiens", "types": "Pathway"}, timeout=2)
        if r.ok:
            for res in r.json().get("results", []):
                paths.add(res["name"])
    except: pass
    try:
        r = requests.get("https://webservice.wikipathways.org/findPathwaysByText",
                         params={"query": gene, "format": "json"}, timeout=2)
        if r.ok:
            for res in r.json().get("result", []):
                if res.get("species") == "Homo sapiens":
                    paths.add(res["name"])
    except: pass
    return list(paths)

def build_pathway_map(gene_evidence):
    print("\n[3] 🚀 Mapping Pathways...")
    pathway_map = defaultdict(set)
    genes_to_map = list(gene_evidence.keys())
    with ThreadPoolExecutor(max_workers=10) as executor:
        results = executor.map(fetch_live_pathways, genes_to_map)
        for gene, paths in zip(genes_to_map, results):
            for p in paths:
                p_clean = p.split(" - ")[0].strip()
                if "Disease" not in p_clean and "metabolism" not in p_clean.lower():
                    pathway_map[p_clean].add(gene)
    print(f"    → Mapped {len(pathway_map)} pathways.")
    return pathway_map

# =========================
# Graph Drawing
# =========================
def draw_graph(compound, gene_evidence, pathway_map, mode_name, outcome, graph_path):
    print("\n[4] 🎨 Generating Graph...")
    G = nx.Graph()
    if pathway_map:
        top_paths = sorted(pathway_map.items(), key=lambda x: len(x[1]), reverse=True)[:20]
        pathway_nodes, gene_nodes = set(), set()
        for p, genes in top_paths:
            pathway_nodes.add(p)
            for g in genes:
                gene_nodes.add(g)
                G.add_edge(g, p)
        plt.figure(figsize=(20, 16))
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_nodes(G, pos, nodelist=list(pathway_nodes), node_color='#7ED957', node_size=2800)
        nx.draw_networkx_nodes(G, pos, nodelist=list(gene_nodes), node_color='#6EC1E4', node_size=1400)
        nx.draw_networkx_edges(G, pos, alpha=0.3, width=1.5)
        labels = {n: n.replace(" ", "\n", 2) if len(n)>15 else n for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight="bold")
        plt.axis('off')
        plt.tight_layout()
    else:
        # fallback simple graph
        G.add_node(compound, color='red')
        top_genes = sorted(gene_evidence.keys(), key=lambda g: len(gene_evidence[g]), reverse=True)[:20]
        for g in top_genes: G.add_edge(compound, g)
        plt.figure(figsize=(14, 12))
        pos = nx.kamada_kawai_layout(G)
        nx.draw_networkx_nodes(G, pos, nodelist=[compound], node_color='#FF6B6B', node_size=3000)
        nx.draw_networkx_nodes(G, pos, nodelist=top_genes, node_color='#6EC1E4', node_size=1500)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        nx.draw_networkx_labels(G, pos, font_size=9, font_weight="bold")
        plt.axis('off')
        plt.tight_layout()
    os.makedirs(os.path.dirname(graph_path), exist_ok=True)
    plt.savefig(graph_path, dpi=300)
    plt.close()
    print(f"    ✔ Graph saved: {graph_path}")

# =========================
# Report + JSON
# =========================
def generate_report(job_id, compound, outcome, mode_name, gene_evidence, pathway_map, paper_ids, results_dir):
    os.makedirs(results_dir, exist_ok=True)
    report_path = os.path.join(results_dir, f"{job_id}.json")
    fingerprint = compute_fingerprint(gene_evidence, pathway_map)
    data = {
        "job_id": job_id,
        "compound": compound,
        "outcome": outcome,
        "mode": mode_name,
        "fingerprint": fingerprint,
        "gene_evidence": gene_evidence,
        "pathway_map": {k: list(v) for k, v in pathway_map.items()},
        "papers": paper_ids
    }
    with open(report_path, "w") as f:
        json.dump(data, f, indent=2)
    print(f"    ✔ JSON saved: {report_path}")
    return report_path

# =========================
# Main job runner
# =========================
def run_relio_job(compound, outcome, mode_name, job_id):
    whitelist = build_gene_whitelist(outcome)
    if not whitelist:
        print("❌ Whitelist failed.")
        return None

    if mode_name.upper() == "FULL":
        gene_evidence, paper_ids = mine_pmc_full(compound, outcome, whitelist, limit=100)
        mode_name_str = "FULL (PMC Full Text)"
    else:
        gene_evidence, paper_ids = mine_abstracts_fast(compound, outcome, whitelist, limit=100)
        mode_name_str = "FAST (Abstracts)"

    if not gene_evidence:
        print("❌ No evidence found.")
        return None

    pathway_map = build_pathway_map(gene_evidence)

    # Paths
    results_dir = os.path.join("python-runner", "results")
    images_dir = os.path.join("python-runner", "images", job_id)
    os.makedirs(images_dir, exist_ok=True)
    graph_path = os.path.join(images_dir, f"{job_id}_graph.png")

    draw_graph(compound, gene_evidence, pathway_map, mode_name_str, outcome, graph_path)
    report_path = generate_report(job_id, compound, outcome, mode_name_str, gene_evidence, pathway_map, paper_ids, results_dir)

    return {"report": report_path, "graph": graph_path}
