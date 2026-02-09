import os
import time
import re
import math
import requests
import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import json

from Bio import Entrez

# =========================
# CONFIGURATION
# =========================
Entrez.email = "your.email@example.com"  # <-- change this
Entrez.tool = "Relio"
NCBI_DELAY = 0.35

# =========================
# UTILITIES
# =========================
def sleep():
    time.sleep(NCBI_DELAY)

def get_context_direction(text, gene):
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts:
        return "neutral"
    idx = starts[0]
    snippet = text[max(0, idx-100): min(len(text), idx+100)].lower()
    if any(w in snippet for w in ["increase", "induce", "upregulat", "activat", "stimulat", "enhance", "promot"]):
        return "up"
    if any(w in snippet for w in ["decrease", "inhibit", "suppress", "reduc", "block", "attenuate", "prevent"]):
        return "down"
    return "neutral"

# =========================
# 1. WHITELIST
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
    except:
        return set()

# =========================
# 2. MINING
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
        if not ids: return {}, []
        print(f"    → Found {len(ids)} abstracts. Downloading...")
        sleep()
        h = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
        text_blob = h.read()
        h.close()
        abstracts = [a.strip() for a in text_blob.split("\n\n") if len(a)>50]
        for i, txt in enumerate(abstracts):
            pid = ids[i] if i < len(ids) else "Unknown"
            for g in genes:
                if re.search(rf"\b{g}\b", txt):
                    direction = get_context_direction(txt, g)
                    gene_evidence[g].append({"id": pid, "direction": direction})
        return gene_evidence, pmids
    except:
        return {}, []

# =========================
# 3. PATHWAY MAPPING
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
    print("\n[3] 🚀 Mapping Pathways (Live API)...")
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
# 4. GRAPH GENERATION
# =========================
def draw_graph(compound, gene_evidence, pathway_map, mode_name, outcome, graph_path):
    print("\n[4] 🎨 Generating Graph...")
    G = nx.Graph()
    if pathway_map:
        top_paths = sorted(pathway_map.items(), key=lambda x: len(x[1]), reverse=True)[:20]
        for p, genes in top_paths:
            for g in genes:
                G.add_edge(g, p)
    else:
        top_genes = sorted(gene_evidence.keys(), key=lambda g: len(gene_evidence[g]), reverse=True)[:20]
        G.add_node(compound)
        for g in top_genes:
            G.add_edge(compound, g)
    plt.figure(figsize=(14,12))
    pos = nx.kamada_kawai_layout(G)
    nx.draw_networkx_nodes(G, pos, node_color='#6EC1E4', node_size=1500)
    nx.draw_networkx_edges(G, pos, alpha=0.4)
    nx.draw_networkx_labels(G, pos, font_size=9, font_weight="bold")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(graph_path, dpi=300)
    plt.close()
    print(f"    ✔ Graph saved: {graph_path}")

# =========================
# 5. REPORT GENERATION
# =========================
def generate_full_report_text(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name):
    lines = []
    lines.append("="*60)
    lines.append(f"Relio REPORT: {compound} + {outcome}")
    lines.append("="*60)
    lines.append(f"Mode: {mode_name}")
    lines.append(f"Documents Analyzed: {len(paper_ids)}")
    lines.append("-"*40)
    lines.append("Gene Evidence Summary:")
    for g, ev in gene_evidence.items():
        hits = len(ev)
        dirs = [e['direction'] for e in ev]
        d_mode = Counter(dirs).most_common(1)[0][0] if dirs else "neutral"
        lines.append(f"{g}: {hits} hits, direction={d_mode}")
    return "\n".join(lines)

# =========================
# 6. MAIN FUNCTION
# =========================
def run_relio_job(compound, outcome, mode_name, job_id):
    whitelist = build_gene_whitelist(outcome)
    if not whitelist:
        print("❌ Whitelist failed.")
        return None

    if mode_name.upper() == "FAST":
        gene_evidence, paper_ids = mine_abstracts_fast(compound, outcome, whitelist)
        mode_name_str = "FAST"
    else:
        # For now only FAST implemented
        gene_evidence, paper_ids = mine_abstracts_fast(compound, outcome, whitelist)
        mode_name_str = "FULL"

    pathway_map = build_pathway_map(gene_evidence)

    # Paths
    results_dir = os.path.join("python-runner", "results")
    os.makedirs(results_dir, exist_ok=True)

    images_dir = os.path.join("python-runner", "images")
    os.makedirs(images_dir, exist_ok=True)

    graph_filename = f"{job_id}.png"
    graph_path = os.path.join(images_dir, graph_filename)

    draw_graph(compound, gene_evidence, pathway_map, mode_name_str, outcome, graph_path)

    report_text = generate_full_report_text(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name_str)

    # Save JSON with actual results
    result_json = {
        "job_id": job_id,
        "compound": compound,
        "outcome": outcome,
        "mode": mode_name_str,
        "report": report_text,
        "graph": graph_filename,
        "gene_evidence": gene_evidence,
        "pathways": {k: list(v) for k,v in pathway_map.items()}
    }

    result_file = os.path.join(results_dir, f"{job_id}.json")
    with open(result_file, "w") as f:
        json.dump(result_json, f, indent=2)

    print(f"\n[5] 📝 Results saved: {result_file}")
    return result_json
