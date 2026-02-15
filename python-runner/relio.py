import time, re, math, os, requests, networkx as nx, matplotlib.pyplot as plt
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from Bio import Entrez

Entrez.email = "relio@app.ai"
NCBI_DELAY = 0.35

def sleep(): time.sleep(NCBI_DELAY)

# ===============================
# CONTEXT DIRECTION
# ===============================
def get_context_direction(text, gene):
    starts = [m.start() for m in re.finditer(rf"\b{gene}\b", text)]
    if not starts: return "neutral"
    idx = starts[0]
    snippet = text[max(0, idx-100):min(len(text), idx+100)].lower()

    if any(w in snippet for w in ["increase","induce","upregulat","activat","stimulat","enhance","promot"]): return "up"
    if any(w in snippet for w in ["decrease","inhibit","suppress","reduc","block","attenuate","prevent"]): return "down"
    return "neutral"

# ===============================
# WHITELIST
# ===============================
def build_gene_whitelist(outcome):
    genes=set()
    try:
        sleep()
        h=Entrez.esearch(db="gene",term=f"{outcome} AND Homo sapiens[Organism]",retmax=200)
        ids=Entrez.read(h)["IdList"]
        h.close()
        if ids:
            sleep()
            h=Entrez.esummary(db="gene",id=",".join(ids))
            data=Entrez.read(h,validate=False)
            h.close()
            for d in data["DocumentSummarySet"]["DocumentSummary"]:
                g=d.get("Name","")
                if g.isupper() and 2<=len(g)<=10: genes.add(g)
    except: pass
    return genes

# ===============================
# ABSTRACT MINING
# ===============================
def mine_abstracts_fast(compound,outcome,genes,limit=100):
    gene_evidence=defaultdict(list)
    try:
        h=Entrez.esearch(db="pubmed",term=f"{compound} AND {outcome}",retmax=limit,sort="relevance")
        ids=Entrez.read(h)["IdList"];h.close()
        if not ids: return {},[]
        sleep()
        h=Entrez.efetch(db="pubmed",id=",".join(ids),rettype="abstract",retmode="text")
        text=h.read();h.close()
        abstracts=[a for a in text.split("\n\n") if len(a)>50]

        for i,txt in enumerate(abstracts):
            pid=ids[i] if i<len(ids) else "NA"
            for g in genes:
                if re.search(rf"\b{g}\b",txt):
                    gene_evidence[g].append({"id":pid,"direction":get_context_direction(txt,g)})
        return gene_evidence,ids
    except:
        return {},[]

# ===============================
# PMC FULL TEXT
# ===============================
def mine_pmc_full(compound,outcome,genes,limit=100):
    gene_evidence=defaultdict(list)
    try:
        h=Entrez.esearch(db="pmc",term=f"{compound} AND {outcome} AND open access[filter]",retmax=limit)
        ids=Entrez.read(h)["IdList"];h.close()
        for i in range(0,len(ids),20):
            sleep()
            h=Entrez.efetch(db="pmc",id=",".join(ids[i:i+20]),retmode="xml")
            arts=Entrez.read(h,validate=False);h.close()
            arts=arts if isinstance(arts,list) else arts.get("pmc-articles",[])
            for art in arts:
                text=str(art)
                for g in genes:
                    if re.search(rf"\b{g}\b",text):
                        gene_evidence[g].append({"id":"PMC","direction":get_context_direction(text,g)})
        return gene_evidence,ids
    except:
        return {},[]

# ===============================
# PATHWAYS
# ===============================
def fetch_live_pathways(gene):
    paths=set()
    try:
        r=requests.get("https://reactome.org/ContentService/search/query",params={"query":gene})
        if r.ok:
            for x in r.json().get("results",[]): paths.add(x["name"])
    except: pass
    return list(paths)

def build_pathway_map(gene_evidence):
    pm=defaultdict(set)
    with ThreadPoolExecutor(max_workers=10) as exe:
        for g,paths in zip(gene_evidence.keys(),exe.map(fetch_live_pathways,gene_evidence.keys())):
            for p in paths: pm[p].add(g)
    return pm

# ===============================
# CONTRADICTIONS
# ===============================
def analyze_contradictions(gene_evidence):
    conflicts=[]
    for g,hits in gene_evidence.items():
        ups=sum(1 for x in hits if x["direction"]=="up")
        downs=sum(1 for x in hits if x["direction"]=="down")
        if ups>0 and downs>0:
            conflicts.append({"gene":g,"up":ups,"down":downs})
    return conflicts

# ===============================
# FINGERPRINT
# ===============================
def compute_fingerprint(gene_evidence,pathway_map):
    up=sum(1 for h in gene_evidence.values() for x in h if x["direction"]=="up")
    down=sum(1 for h in gene_evidence.values() for x in h if x["direction"]=="down")
    if up>down*1.5: d="Activation"
    elif down>up*1.5: d="Inhibition"
    else: d="Balanced"
    return {"direction":d,"entropy":len(pathway_map)}

# ===============================
# GRAPH
# ===============================
def draw_graph(pathway_map,job_id):
    if not pathway_map: return None
    G=nx.Graph()
    for p,gs in pathway_map.items():
        for g in gs: G.add_edge(g,p)
    plt.figure(figsize=(16,14))
    nx.draw(G,with_labels=False,node_size=500)
    os.makedirs("python-runner/images",exist_ok=True)
    path=f"python-runner/images/{job_id}_graph.png"
    plt.savefig(path,dpi=200)
    plt.close()
    return path

# ===============================
# MAIN
# ===============================
def run_relio(compound,outcome,mode,job_id):
    whitelist=build_gene_whitelist(outcome)
    if not whitelist: return {"error":"whitelist"}

    if mode=="FAST":
        ev,docs=mine_abstracts_fast(compound,outcome,whitelist)
    else:
        ev,docs=mine_pmc_full(compound,outcome,whitelist)

    if not ev: return {"error":"no evidence"}

    pmap=build_pathway_map(ev)
    fp=compute_fingerprint(ev,pmap)
    conflicts=analyze_contradictions(ev) if mode=="FULL" else []
    img=draw_graph(pmap,job_id)

    return {
        "compound":compound,
        "outcome":outcome,
        "mode":mode,
        "fingerprint":fp,
        "genes":ev,
        "pathways":{k:list(v) for k,v in pmap.items()},
        "contradictions":conflicts,
        "documents":docs,
        "images":{"graph":os.path.basename(img) if img else None}
    }
