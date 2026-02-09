import time, re, math, os, requests
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import networkx as nx
import matplotlib.pyplot as plt

# =========================
# 1. CONFIG
# =========================
try:
    from Bio import Entrez
except ImportError:
    print("Install biopython: pip install biopython networkx matplotlib requests")
    exit()

Entrez.email = "your.email@example.com"
Entrez.tool = "Relio"
NCBI_DELAY = 0.35

# =========================
# 2. UTILITIES
# =========================
def sleep(): time.sleep(NCBI_DELAY)

def get_context_direction(text, gene):
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts: return "neutral"
    idx = starts[0]
    snippet = text[max(0, idx-100): idx+100].lower()
    if any(w in snippet for w in ["increase","induce","upregulat","activat","stimulat","enhance"]): return "up"
    if any(w in snippet for w in ["decrease","inhibit","suppress","reduc","block","attenuate"]): return "down"
    return "neutral"

# =========================
# 3. WHITELIST
# =========================
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
                g = d.get("Name","")
                if g.isupper() and 2<=len(g)<=10: genes.add(g)
        print(f"    → Whitelist: {len(genes)} genes.")
        return genes
    except: return set()

# =========================
# 4. MINING
# =========================
def mine_abstracts_fast(compound, outcome, genes, limit=100):
    print(f"\n[2] ⚡ FAST MODE: Fetching up to {limit} PubMed Abstracts...")
    gene_evidence = defaultdict(list)
    try:
        h = Entrez.esearch(db="pubmed", term=f"{compound} AND {outcome}", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        if not ids: return {}, []
        sleep()
        h = Entrez.efetch(db="pubmed", id=",".join(ids), rettype="abstract", retmode="text")
        text_blob = h.read()
        h.close()
        abstracts = [a.strip() for a in text_blob.split("\n\n") if len(a)>50]
        for i, txt in enumerate(abstracts):
            pid = ids[i] if i<len(ids) else "Unknown"
            for g in genes:
                if re.search(rf"\b{g}\b", txt):
                    direction = get_context_direction(txt, g)
                    gene_evidence[g].append({"id": pid, "direction": direction})
        return gene_evidence, ids
    except: return {}, []

def mine_pmc_full(compound, outcome, genes, limit=100):
    print(f"\n[2] 🤿 FULL MODE: Searching up to {limit} PMC Full Text Articles...")
    gene_evidence = defaultdict(list)
    try:
        h = Entrez.esearch(db="pmc", term=f"{compound} AND {outcome} AND open access[filter]", retmax=limit)
        ids = Entrez.read(h)["IdList"]
        h.close()
        if not ids: return {}, []
        batch_size = 20
        for i in range(0,len(ids),batch_size):
            sleep()
            h = Entrez.efetch(db="pmc", id=",".join(ids[i:i+batch_size]), retmode="xml")
            articles = Entrez.read(h, validate=False)
            h.close()
            art_list = articles if isinstance(articles,list) else articles.get("pmc-articles",[])
            if isinstance(art_list, dict): art_list=[art_list]
            for art in art_list:
                current_pmc="Unknown"
                try:
                    for aid in art.get('front',{}).get('article-meta',{}).get('article-id',[]):
                        if aid.attributes.get('pub-id-type')=="pmc": current_pmc="PMC"+str(aid)
                except: pass
                full_text = str(art)
                for g in genes:
                    if re.search(rf"\b{g}\b", full_text): 
                        direction = get_context_direction(full_text,g)
                        gene_evidence[g].append({"id": current_pmc,"direction":direction})
        return gene_evidence, ids
    except: return {}, []

# =========================
# 5. PATHWAY MAPPING
# =========================
def fetch_live_pathways(gene):
    paths = set()
    try:
        r = requests.get("https://reactome.org/ContentService/search/query", params={"query":gene,"species":"Homo sapiens","types":"Pathway"}, timeout=2)
        if r.ok:
            for res in r.json().get("results",[]): paths.add(res["name"])
    except: pass
    return list(paths)

def build_pathway_map(gene_evidence):
    print("\n[3] 🚀 Mapping Pathways...")
    pathway_map = defaultdict(set)
    genes_to_map = list(gene_evidence.keys())
    with ThreadPoolExecutor(max_workers=10) as executor:
        results = executor.map(fetch_live_pathways, genes_to_map)
        for gene, paths in zip(genes_to_map, results):
            for p in paths: pathway_map[p].add(gene)
    print(f"    → Mapped {len(pathway_map)} pathways.")
    return pathway_map

# =========================
# 6. GRAPH
# =========================
def draw_graph(compound, gene_evidence, pathway_map, mode_name, outcome, image_dir):
    print("\n[4] 🎨 Generating Graph...")
    G = nx.Graph()
    for p, genes in pathway_map.items():
        for g in genes: G.add_edge(g,p)
    plt.figure(figsize=(12,10))
    pos = nx.kamada_kawai_layout(G)
    nx.draw(G,pos, with_labels=True, node_color="#6EC1E4")
    plt.title(f"{compound} + {outcome} ({mode_name})")
    plt.axis("off")
    os.makedirs(image_dir, exist_ok=True)
    plt.savefig(os.path.join(image_dir,"graph.png"), dpi=300)
    print(f"    ✔ Graph saved: {os.path.join(image_dir,'graph.png')}")

# =========================
# 7. REPORT
# =========================
def generate_full_report(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name, result_dir):
    print("\n[5] 📝 Generating Report...")
    lines = [f"Relio Report: {compound} + {outcome} ({mode_name})","="*60]
    lines.append("Gene Evidence:")
    for g, ev in gene_evidence.items():
        hits = len(ev)
        dirs = [e['direction'] for e in ev]
        d_mode = Counter(dirs).most_common(1)[0][0]
        lines.append(f"{g}: {hits} hits ({d_mode})")
    os.makedirs(result_dir, exist_ok=True)
    report_path = os.path.join(result_dir,"report.txt")
    with open(report_path,"w") as f:
        f.write("\n".join(lines))
    print(f"    ✔ Report saved: {report_path}")

# =========================
# 8. MAIN FUNCTION
# =========================
def run_relio_job(data, job_id, result_dir, image_dir):
    compound = data["compound"]
    outcome = data["outcome"]
    mode = data["mode"].upper()
    whitelist = build_gene_whitelist(outcome)
    if not whitelist:
        print("❌ Whitelist failed")
        return
    if mode=="FAST":
        ev,paper_ids = mine_abstracts_fast(compound,outcome,whitelist)
    else:
        ev,paper_ids = mine_pmc_full(compound,outcome,whitelist)
    if not ev:
        print("❌ No evidence found")
        return
    pmap = build_pathway_map(ev)
    draw_graph(compound, ev, pmap, mode, outcome, image_dir)
    generate_full_report(compound, outcome, ev, pmap, paper_ids, mode, result_dir)
