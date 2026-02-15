import time, re, os, requests
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez

Entrez.email = "relio@github"
NCBI_DELAY = 0.35

def sleep():
    time.sleep(NCBI_DELAY)

def get_context_direction(text, gene):
    t = text.lower()
    if "increase" in t or "activate" in t: return "up"
    if "decrease" in t or "inhibit" in t: return "down"
    return "neutral"

def analyze_contradictions(genes):
    out = []
    for g, h in genes.items():
        up = sum(1 for x in h if x["direction"] == "up")
        down = sum(1 for x in h if x["direction"] == "down")
        if up and down:
            out.append({"gene": g, "up": up, "down": down})
    return out

def compute_fingerprint(genes, pathways):
    up = sum(1 for h in genes.values() for x in h if x["direction"] == "up")
    down = sum(1 for h in genes.values() for x in h if x["direction"] == "down")
    return {
        "direction": "Activation" if up >= down else "Inhibition",
        "entropy": len(pathways)
    }

def build_gene_whitelist(outcome):
    sleep()
    h = Entrez.esearch(db="gene", term=f"{outcome} AND Homo sapiens", retmax=50)
    ids = Entrez.read(h)["IdList"]
    h.close()
    return set(ids)

def mine_abstracts_fast(compound, outcome, genes):
    sleep()
    h = Entrez.esearch(db="pubmed", term=f"{compound} AND {outcome}", retmax=20)
    ids = Entrez.read(h)["IdList"]
    h.close()
    if not ids:
        return {}, []

    ev = defaultdict(list)
    for pid in ids:
        ev[pid].append({"direction": "up"})
    return ev, ids

def mine_pmc_full(*args):
    return mine_abstracts_fast(*args)

def build_pathway_map(genes):
    p = defaultdict(set)
    for g in genes:
        p["Immune"].add(g)
    return p

def draw_graph(compound, genes, pathways, mode, outcome):
    BASE = os.path.dirname(__file__)
    IMG_DIR = os.path.join(BASE, "images")
    os.makedirs(IMG_DIR, exist_ok=True)

    G = nx.Graph()
    for p, gs in pathways.items():
        for g in gs:
            G.add_edge(g, p)

    if not G.nodes:
        G.add_node(compound)

    plt.figure(figsize=(12, 10))
    pos = nx.spring_layout(G, k=0.8)
    nx.draw(G, pos, with_labels=True, node_size=1200, font_size=8)

    outfile = os.path.join(IMG_DIR, "litmap_maze.png")
    plt.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close()

    print("GRAPH GENERATED:", outfile)

def generate_full_report(compound, outcome, genes, pathways, docs, mode):
    BASE = os.path.dirname(__file__)
    RES_DIR = os.path.join(BASE, "results")
    os.makedirs(RES_DIR, exist_ok=True)

    with open(os.path.join(RES_DIR, "litmap_report.txt"), "w") as f:
        f.write(f"{compound} vs {outcome}\nDocuments: {len(docs)}\n")
