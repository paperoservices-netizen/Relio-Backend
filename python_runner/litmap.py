import time, re, os, requests, math
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez

Entrez.email = "relio@github"
NCBI_DELAY = 0.35

def sleep():
    time.sleep(NCBI_DELAY)

def get_context_direction(text, gene):
    """Scans for up/down regulation keywords."""
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts: return "neutral"
    idx = starts[0]
    snippet = text[max(0, idx-100): min(len(text), idx+100)].lower()

    if any(w in snippet for w in ["increase", "induce", "upregulat", "activat", "stimulat", "enhance", "promot"]): return "up"
    if any(w in snippet for w in ["decrease", "inhibit", "suppress", "reduc", "block", "attenuate", "prevent"]): return "down"
    return "neutral"

def analyze_contradictions(gene_evidence):
    """Identifies genes with conflicting evidence."""
    conflicts = []
    for gene, hits in gene_evidence.items():
        dirs = [h['direction'] for h in hits]
        up = dirs.count('up')
        down = dirs.count('down')
        total = up + down
        if total > 1:
            minority = min(up, down)
            if minority > 0 and (minority / total) > 0.2:
                conflicts.append({
                    "gene": gene, "up": up, "down": down, "consensus": "Contradictory"
                })
    return conflicts

def compute_fingerprint(gene_evidence, pathway_map):
    up = sum(1 for hits in gene_evidence.values() for h in hits if h['direction']=='up')
    down = sum(1 for hits in gene_evidence.values() for h in hits if h['direction']=='down')

    if up > down * 1.5: direction = "Predominant Activation"
    elif down > up * 1.5: direction = "Predominant Inhibition"
    else: direction = "Balanced Modulation"

    path_counts = [len(genes) for genes in pathway_map.values()]
    total_conn = sum(path_counts)
    entropy = 0.0
    if total_conn > 0:
        probs = [c / total_conn for c in path_counts]
        entropy = -sum(p * math.log2(p) for p in probs)

    return {"direction": direction, "entropy": entropy}

def build_gene_whitelist(outcome):
    print(f"\n[1] 🧠 Building Gene Whitelist for '{outcome}'...")
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
                if g.isupper() and 2 <= len(g) <= 10: genes.add(g)
        print(f"    → Whitelist: {len(genes)} verified genes.")
        return genes
    except: 
        print("    → Whitelist failed.")
        return set()

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
    except: return {}, []

def mine_pmc_full(compound, outcome, genes, limit=100):
    print(f"\n[2] 🤿 FULL MODE: Searching up to {limit} PMC Full Text Articles...")
    gene_evidence = defaultdict(list)
    pmc_list = []
    try:
        h = Entrez.esearch(db="pmc", term=f"{compound} AND {outcome} AND open access[filter]", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        pmc_list = ids
        if not ids: return {}, []
        print(f"    → Found {len(ids)} full-text XMLs. Downloading...")

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
                            if aid.attributes.get('pub-id-type') == 'pmc': current_pmc = "PMC" + str(aid)
                    except: pass
                    full_text = str(art)
                    for g in genes:
                        if re.search(rf"[^a-zA-Z]{g}[^a-zA-Z]", full_text):
                            direction = get_context_direction(full_text, g)
                            gene_evidence[g].append({"id": current_pmc, "direction": direction})
            except: continue
        return gene_evidence, pmc_list
    except: return {}, []

def fetch_live_pathways(gene):
    paths = set()
    try:
        r = requests.get("https://reactome.org/ContentService/search/query",
                       params={"query": gene, "species": "Homo sapiens", "types": "Pathway"}, timeout=2)
        if r.ok:
            for res in r.json().get("results", []): paths.add(res["name"])
    except: pass
    try:
        r = requests.get("https://webservice.wikipathways.org/findPathwaysByText",
                       params={"query": gene, "format": "json"}, timeout=2)
        if r.ok:
            for res in r.json().get("result", []):
                if res.get("species") == "Homo sapiens": paths.add(res["name"])
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

def draw_graph(compound, gene_evidence, pathway_map, mode_name, outcome):
    if "FAST" in mode_name:
        print("\n[4] ⏩ Graph generation skipped (Fast Mode).")
        return

    print("\n[4] 🎨 Generating Maze-Like Graph...")
    G = nx.Graph()
    if pathway_map:
        top_paths = sorted(pathway_map.items(), key=lambda x: len(x[1]), reverse=True)[:20]
        pathway_nodes = set()
        gene_nodes = set()
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

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='Target Gene', markerfacecolor='#6EC1E4', markersize=15),
            Line2D([0], [0], marker='o', color='w', label='Biological Pathway', markerfacecolor='#7ED957', markersize=15)
        ]
        plt.legend(handles=legend_elements, loc='upper left', fontsize=12, frameon=True)
        plt.suptitle(f"LitMap Analysis: {compound} + {outcome}", fontsize=16, fontweight='bold', y=0.95)

    else:
        print("    ⚠️ Using Gene-Star Graph (No pathways found)")
        top_genes = sorted(gene_evidence.keys(), key=lambda g: len(gene_evidence[g]), reverse=True)[:20]
        G.add_node(compound, color='red')
        for g in top_genes: G.add_edge(compound, g)

        plt.figure(figsize=(14, 12))
        pos = nx.kamada_kawai_layout(G)

        nx.draw_networkx_nodes(G, pos, nodelist=[compound], node_color='#FF6B6B', node_size=3000)
        nx.draw_networkx_nodes(G, pos, nodelist=top_genes, node_color='#6EC1E4', node_size=1500)
        nx.draw_networkx_edges(G, pos, alpha=0.4)
        nx.draw_networkx_labels(G, pos, font_size=9, font_weight="bold")

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label='Compound', markerfacecolor='#FF6B6B', markersize=15),
            Line2D([0], [0], marker='o', color='w', label='Gene Target', markerfacecolor='#6EC1E4', markersize=15)
        ]
        plt.legend(handles=legend_elements, loc='lower right', fontsize=12)

    plt.axis('off')
    plt.tight_layout()
    plt.savefig("litmap_maze.png", dpi=300)
    plt.close()
    print("    ✔ Graph saved as 'litmap_maze.png'")

def generate_full_report(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name):
    print("\n[5] 📝 Generating Report...")
    fingerprint = compute_fingerprint(gene_evidence, pathway_map)

    conflicts = []
    if "FULL" in mode_name:
        conflicts = analyze_contradictions(gene_evidence)

    gene_to_paths = defaultdict(list)
    for p, genes in pathway_map.items():
        for g in genes: gene_to_paths[g].append(p)

    lines = []
    lines.append("="*80)
    lines.append(f"LITMAP V0.1 REPORT: {compound.upper()} + {outcome.upper()}")
    lines.append("="*80)
    lines.append(f"Mode              : {mode_name}")
    lines.append(f"Docs Analyzed     : {len(paper_ids)} (Strict Limit: 100)")
    lines.append(f"Net Direction     : {fingerprint['direction']}")
    lines.append("-" * 80 + "\n")

    if "FULL" in mode_name:
        lines.append("1. CONTRADICTION & CONTEXT ANALYSIS")
        lines.append("-" * 40)
        if conflicts:
            lines.append("Significant conflicting evidence found for the following genes:")
            for c in conflicts:
                lines.append(f"- {c['gene']}: {c['up']} UP studies vs {c['down']} DOWN studies.")
        else:
            lines.append("No significant contradictions found in the full text analysis.\n")

    lines.append("2. GENE ASSOCIATION LIST (Complete Pathway Mapping)")
    lines.append("-" * 80)
    lines.append(f"{'GENE':<8} | {'HITS':<4} | {'DIR':<7} | {'FULL PATHWAY LIST'}")
    lines.append("-" * 80)

    top_genes_all = sorted(gene_evidence.items(), key=lambda x: len(x[1]), reverse=True)
    for g, ev in top_genes_all:
        hits = len(ev)
        dirs = [e['direction'] for e in ev]
        d_mode = Counter(dirs).most_common(1)[0][0]
        if any(c['gene'] == g for c in conflicts): d_mode = "MIXED"

        my_paths = gene_to_paths.get(g, [])
        path_str = ", ".join(my_paths) if my_paths else "No specific pathways mapped"
        lines.append(f"{g:<8} | {hits:<4} | {d_mode:<7} | {path_str}")
    lines.append("-" * 80 + "\n")

    lines.append("3. SOURCE DOCUMENT LIST (Clickable Links)")
    lines.append("-" * 40)
    unique_ids = sorted(list(set(paper_ids)))
    for pid in unique_ids:
        pid_str = str(pid)
        if mode_name.startswith("FULL"):
            if "PMC" not in pid_str and pid_str != "Unknown": pid_str = "PMC" + pid_str
            if "PMC" in pid_str:
                lines.append(f"- https://www.ncbi.nlm.nih.gov/pmc/articles/{pid_str}/")
        else:
            lines.append(f"- https://pubmed.ncbi.nlm.nih.gov/{pid_str}/")

    final_report_str = "\n".join(lines)

    with open("litmap_report.txt", "w") as f:
        f.write(final_report_str)
    print("    ✔ Report saved as 'litmap_report.txt'")

    print("\n" + "="*30 + " TERMINAL REPORT OUTPUT " + "="*30 + "\n")
    print(final_report_str)
    print("\n" + "="*80 + "\n")
