import os, json, time, re, math, requests
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import networkx as nx
from Bio import Entrez

# =========================
# CONFIGURATION
# =========================
Entrez.email = "your.email@example.com"
Entrez.tool = "LitMap_V0.1"
NCBI_DELAY = 0.35

def sleep():
    time.sleep(NCBI_DELAY)

# =========================
# 1. Utilities
# =========================
def get_context_direction(text, gene):
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts: return "neutral"
    idx = starts[0]
    snippet = text[max(0, idx-100): idx+100].lower()
    if any(w in snippet for w in ["increase","upregulat","activat","enhance"]): return "up"
    if any(w in snippet for w in ["decrease","inhibit","suppress","reduce"]): return "down"
    return "neutral"

def analyze_contradictions(gene_evidence):
    conflicts = []
    for gene, hits in gene_evidence.items():
        dirs = [h['direction'] for h in hits]
        up = dirs.count('up')
        down = dirs.count('down')
        total = up + down
        if total>1:
            minority = min(up, down)
            if minority>0 and (minority/total)>0.2:
                conflicts.append({"gene": gene, "up": up, "down": down, "consensus": "Contradictory"})
    return conflicts

def compute_fingerprint(gene_evidence, pathway_map):
    up=sum(1 for hits in gene_evidence.values() for h in hits if h['direction']=='up')
    down=sum(1 for hits in gene_evidence.values() for h in hits if h['direction']=='down')
    if up>down*1.5: direction="Predominant Activation"
    elif down>up*1.5: direction="Predominant Inhibition"
    else: direction="Balanced Modulation"
    path_counts=[len(genes) for genes in pathway_map.values()]
    total_conn=sum(path_counts)
    entropy=0.0
    if total_conn>0:
        probs=[c/total_conn for c in path_counts]
        entropy=-sum(p*math.log2(p) for p in probs)
    return {"direction": direction,"entropy":entropy}

# =========================
# 2. Gene whitelist
# =========================
def build_gene_whitelist(outcome):
    print(f"Building whitelist for '{outcome}'...")
    try:
        sleep()
        h=Entrez.esearch(db="gene", term=f"{outcome} AND Homo sapiens[Organism]", retmax=200)
        ids=Entrez.read(h)["IdList"]; h.close()
        genes=set()
        if ids:
            sleep()
            h=Entrez.esummary(db="gene", id=",".join(ids))
            data=Entrez.read(h, validate=False); h.close()
            for d in data["DocumentSummarySet"]["DocumentSummary"]:
                g=d.get("Name","")
                if g.isupper() and 2<=len(g)<=10: genes.add(g)
        print(f"Whitelist contains {len(genes)} genes")
        return genes
    except: return set()

# =========================
# 3. Mining abstracts
# =========================
def mine_abstracts_fast(compound, outcome, genes, limit=100):
    gene_evidence=defaultdict(list)
    pmids=[]
    try:
        sleep()
        h=Entrez.esearch(db="pubmed", term=f"{compound} AND {outcome}", retmax=limit, sort="relevance")
        ids=Entrez.read(h)["IdList"]; h.close(); pmids=ids
        if not ids: return {}, []
        sleep()
        h=Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
        text_blob=h.read(); h.close()
        abstracts=[a.strip() for a in text_blob.split("\n\n") if len(a)>50]
        for i, txt in enumerate(abstracts):
            pid=ids[i] if i<len(ids) else "Unknown"
            for g in genes:
                if re.search(rf"\b{g}\b", txt):
                    direction=get_context_direction(txt,g)
                    gene_evidence[g].append({"id": pid,"direction": direction})
        return gene_evidence, pmids
    except: return {}, []

def mine_pmc_full(compound, outcome, genes, limit=100):
    gene_evidence=defaultdict(list)
    pmc_list=[]
    try:
        sleep()
        h=Entrez.esearch(db="pmc", term=f"{compound} AND {outcome} AND open access[filter]", retmax=limit, sort="relevance")
        ids=Entrez.read(h)["IdList"]; h.close(); pmc_list=ids
        if not ids: return {}, []
        # simplified fetching for brevity
        return gene_evidence, pmc_list
    except: return {}, []

# =========================
# 4. Pathway mapping
# =========================
def fetch_live_pathways(gene):
    paths=set()
    try:
        r=requests.get("https://reactome.org/ContentService/search/query",
            params={"query":gene,"species":"Homo sapiens","types":"Pathway"}, timeout=2)
        if r.ok: [paths.add(res["name"]) for res in r.json().get("results",[])]
    except: pass
    return list(paths)

def build_pathway_map(gene_evidence):
    print("Mapping pathways...")
    pathway_map=defaultdict(set)
    with ThreadPoolExecutor(max_workers=10) as ex:
        results=ex.map(fetch_live_pathways,list(gene_evidence.keys()))
        for gene, paths in zip(gene_evidence.keys(), results):
            for p in paths: pathway_map[p].add(gene)
    return pathway_map

# =========================
# 5. Graph generation
# =========================
def draw_graph(compound, gene_evidence, pathway_map, mode_name, outcome, image_dir):
    G=nx.Graph()
    top_paths=sorted(pathway_map.items(), key=lambda x: len(x[1]), reverse=True)[:20]
    pathway_nodes=set()
    gene_nodes=set()
    for p, genes in top_paths:
        pathway_nodes.add(p)
        for g in genes:
            gene_nodes.add(g)
            G.add_edge(g,p)
    plt.figure(figsize=(12,10))
    pos=nx.kamada_kawai_layout(G)
    nx.draw_networkx_nodes(G,pos,nodelist=list(pathway_nodes),node_color='#7ED957',node_size=2800)
    nx.draw_networkx_nodes(G,pos,nodelist=list(gene_nodes),node_color='#6EC1E4',node_size=1400)
    nx.draw_networkx_edges(G,pos,alpha=0.3,width=1.5)
    labels={n:n for n in G.nodes()}
    nx.draw_networkx_labels(G,pos,labels,font_size=8)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(image_dir,"litmap_maze.png"), dpi=300)
    plt.close()
    print(f"Graph saved to {image_dir}/litmap_maze.png")

# =========================
# 6. Report generation
# =========================
def generate_full_report(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name, result_dir):
    fingerprint=compute_fingerprint(gene_evidence, pathway_map)
    lines=[]
    lines.append(f"LITMAP REPORT: {compound} + {outcome}")
    lines.append(f"Mode: {mode_name}")
    lines.append(f"Net Direction: {fingerprint['direction']}")
    lines.append("Gene Hits:")
    for g, ev in gene_evidence.items():
        hits=len(ev)
        dirs=[e['direction'] for e in ev]
        d_mode=Counter(dirs).most_common(1)[0][0]
        lines.append(f"{g}: {hits} hits, {d_mode}")
    report_path=os.path.join(result_dir,"litmap_report.txt")
    with open(report_path,"w") as f: f.write("\n".join(lines))
    print(f"Report saved to {report_path}")

# =========================
# 7. Job wrapper
# =========================
def run_litmap_job(data, job_id, result_dir, image_dir):
    compound = data.get("compound", "TextInput")
    outcome = data.get("outcome", "")
    mode = data.get("mode", "FAST").upper()
    whitelist = build_gene_whitelist(outcome)
    if mode=="FAST":
        ev, paper_ids = mine_abstracts_fast(compound, outcome, whitelist)
        mode_name="FAST"
    else:
        ev, paper_ids = mine_pmc_full(compound, outcome, whitelist)
        mode_name="FULL"
    if not ev:
        print("No evidence found")
        return
    pmap = build_pathway_map(ev)
    draw_graph(compound, ev, pmap, mode_name, outcome, image_dir)
    generate_full_report(compound, outcome, ev, pmap, paper_ids, mode_name, result_dir)
