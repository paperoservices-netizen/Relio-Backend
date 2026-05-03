"""
RELIO V0.1 - GITHUB ACTIONS RUNNER
==================================
Executes the same analysis as Collab code
Outputs JSON results for Cloudflare Worker consumption
"""

import time
import re
import math
import json
import os
import sys
import requests
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
import lxml.etree as ET
from tqdm import tqdm

try:
    from Bio import Entrez
except ImportError:
    print("CRITICAL: Install biopython: pip install biopython networkx matplotlib requests")
    sys.exit(1)

# =========================
# CONFIGURATION
# =========================
Entrez.email = "robio.ra.bt@gmail.com"
Entrez.api_key = "0968ae56e9a676e026f2fd87dcc17a9f8009"
Entrez.tool = "Relio_V0.1_GithubActions"
NCBI_DELAY = 0.11

# Job ID and Paths
JOB_ID = os.environ.get("JOB_ID", "unknown")
RESULTS_DIR = os.path.abspath("./results")
IMAGES_DIR = os.path.abspath("./images")
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(IMAGES_DIR, exist_ok=True)

# =========================
# UTILITIES
# =========================
def sleep_safe():
    """Respect NCBI rate limits"""
    time.sleep(NCBI_DELAY)

def get_context_direction(text, gene):
    """Scans for up/down regulation keywords around the gene."""
    starts = [m.start() for m in re.finditer(rf"\b{re.escape(gene)}\b", text, re.IGNORECASE)]
    if not starts:
        return "neutral"
    
    idx = starts[0]
    snippet = text[max(0, idx-150): min(len(text), idx+150)].lower()
    
    up_words = ["increase", "induce", "upregulat", "activat", "stimulat", "enhance", "promot", "up-regulat", "elevat"]
    down_words = ["decrease", "inhibit", "suppress", "reduc", "block", "attenuate", "prevent", "down-regulat", "knockdown", "silenc"]
    
    if any(w in snippet for w in up_words):
        return "up"
    if any(w in snippet for w in down_words):
        return "down"
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
                    "gene": gene,
                    "up": up,
                    "down": down,
                    "consensus": "Contradictory"
                })
    return conflicts

def compute_fingerprint(gene_evidence, pathway_map):
    """Compute overall direction and entropy of the pathway network."""
    up = sum(1 for hits in gene_evidence.values() for h in hits if h['direction'] == 'up')
    down = sum(1 for hits in gene_evidence.values() for h in hits if h['direction'] == 'down')
    
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
        entropy = -sum(p * math.log2(p) for p in probs if p > 0)
    
    return {"direction": direction, "entropy": round(entropy, 4)}

# =========================
# GENE WHITELIST
# =========================
def build_gene_whitelist(outcome):
    """Build list of verified human genes related to outcome."""
    print(f"\n[1] 🧠 Building Gene Whitelist for '{outcome}'...")
    try:
        sleep_safe()
        h = Entrez.esearch(db="gene", term=f"{outcome} AND Homo sapiens[Organism]", retmax=200)
        ids = Entrez.read(h)["IdList"]
        h.close()
        
        genes = set()
        if ids:
            sleep_safe()
            h = Entrez.esummary(db="gene", id=",".join(ids))
            data = Entrez.read(h, validate=False)
            h.close()
            for d in data["DocumentSummarySet"]["DocumentSummary"]:
                g = d.get("Name", "")
                if g.isupper() and 2 <= len(g) <= 10:
                    genes.add(g)
        
        print(f"    → Whitelist: {len(genes)} verified genes.")
        return genes
    except Exception as e:
        print(f"    ⚠️ Whitelist error: {e}")
        return set()

# =========================
# MINING ENGINES
# =========================
def mine_abstracts_fast(compound, outcome, genes, limit=100):
    """FAST MODE: Mine PubMed abstracts (max 100)"""
    print(f"\n[2] ⚡ FAST MODE: Fetching up to {limit} PubMed Abstracts...")
    gene_evidence = defaultdict(list)
    paper_ids = []
    
    try:
        h = Entrez.esearch(db="pubmed", term=f"{compound} AND {outcome}", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        
        if not ids:
            print("    → No abstracts found.")
            return {}, []
        
        paper_ids = ids
        print(f"    → Found {len(ids)} abstracts. Downloading...")
        
        sleep_safe()
        h = Entrez.efetch(db="pubmed", id=",".join(ids), retmode="xml")
        records = Entrez.read(h)
        h.close()
        
        for article in records.get('PubmedArticle', []):
            try:
                pid = str(article['MedlineCitation']['PMID'])
                abstract_parts = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                full_abstract = " ".join([str(part) for part in abstract_parts])
                
                if not full_abstract:
                    continue
                
                for g in genes:
                    if re.search(rf"\b{re.escape(g)}\b", full_abstract, re.IGNORECASE):
                        direction = get_context_direction(full_abstract, g)
                        gene_evidence[g].append({"id": pid, "direction": direction})
            except KeyError:
                continue
        
        return dict(gene_evidence), paper_ids
    except Exception as e:
        print(f"    ❌ Error in fast mode: {e}")
        return {}, []

def mine_pmc_full(compound, outcome, genes, limit=100):
    """FULL MODE: Mine PMC full-text articles (max 100)"""
    print(f"\n[2] 🤿 FULL MODE: Searching up to {limit} PMC Full Text Articles...")
    gene_evidence = defaultdict(list)
    pmc_list = []
    
    try:
        h = Entrez.esearch(db="pmc", term=f"{compound} AND {outcome} AND open access[filter]", retmax=limit, sort="relevance")
        ids = Entrez.read(h)["IdList"]
        h.close()
        pmc_list = ids
        
        if not ids:
            print("    → No full-text articles found.")
            return {}, []
        
        print(f"    → Found {len(ids)} full-text XMLs. Downloading in safe batches...")
        
        batch_size = 5
        articles_processed = 0
        
        for i in tqdm(range(0, len(ids), batch_size), desc="Mining PMC XMLs", unit="batch", colour="green"):
            batch_ids = ids[i:i+batch_size]
            try:
                sleep_safe()
                h = Entrez.efetch(db="pmc", id=",".join(batch_ids), retmode="xml")
                xml_data = h.read()
                h.close()
                
                if isinstance(xml_data, bytes):
                    xml_data = xml_data.decode('utf-8', errors='ignore')
                
                root = ET.fromstring(xml_data)
                
                if root.tag == 'ERROR' or root.find('.//ERROR') is not None:
                    err_msg = "".join(root.itertext())
                    print(f"    ⚠️ NCBI API Warning: {err_msg.strip()}")
                    continue
                
                for article in root.findall('.//article'):
                    articles_processed += 1
                    current_pmc = "Unknown"
                    
                    for article_id in article.findall('.//article-id'):
                        if article_id.get('pub-id-type') == 'pmc':
                            current_pmc = "PMC" + (article_id.text or "")
                            break
                    
                    text_chunks = []
                    for section in article.findall('.//abstract') + article.findall('.//body'):
                        text_chunks.append(" ".join(section.itertext()))
                    
                    clean_full_text = " ".join(text_chunks)
                    
                    if not clean_full_text.strip():
                        continue
                    
                    for g in genes:
                        if re.search(rf"\b{re.escape(g)}\b", clean_full_text, re.IGNORECASE):
                            direction = get_context_direction(clean_full_text, g)
                            gene_evidence[g].append({"id": current_pmc, "direction": direction})
                
                print(f"    ⏳ Parsed {articles_processed}/{len(ids)} articles...")
                        
            except ET.ParseError as e:
                print(f"    ⚠️ XML Parse Error: {e}")
                continue
            except Exception as e:
                print(f"    ⚠️ Batch error: {e}")
                continue
        
        return dict(gene_evidence), pmc_list
        
    except Exception as e:
        print(f"    ❌ Fatal error in full mode: {e}")
        return {}, []

# =========================
# PATHWAY MAPPING
# =========================
def fetch_live_pathways(gene):
    """Fetch pathways for a gene from Reactome and WikiPathways APIs"""
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
    except:
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
    except:
        pass
    
    return list(paths)

def build_pathway_map(gene_evidence):
    """Build pathway network from live API calls"""
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
    return dict(pathway_map)

# =========================
# GRAPH GENERATION
# =========================
def draw_graph(compound, gene_evidence, pathway_map, mode_name, outcome):
    """Generate and save maze-like graph with legend"""
    if "FAST" in mode_name:
        print("\n[4] ⏩ Graph generation skipped (Fast Mode).")
        return None
    
    print("\n[4] 🎨 Generating Maze-Like Graph...")
    graph_filename = None
    
    try:
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
            
            labels = {n: n.replace(" ", "\n", 2) if len(n) > 15 else n for n in G.nodes()}
            nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight="bold")
            
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', label='Target Gene', markerfacecolor='#6EC1E4', markersize=15),
                Line2D([0], [0], marker='o', color='w', label='Biological Pathway', markerfacecolor='#7ED957', markersize=15)
            ]
            plt.legend(handles=legend_elements, loc='upper left', fontsize=12, frameon=True)
            plt.suptitle(f"RELIO Analysis: {compound} + {outcome}", fontsize=16, fontweight='bold', y=0.95)
        
        else:
            print("    ⚠️ Using Gene-Star Graph (No pathways found)")
            top_genes = sorted(gene_evidence.keys(), key=lambda g: len(gene_evidence[g]), reverse=True)[:20]
            G.add_node(compound)
            for g in top_genes:
                G.add_edge(compound, g)
            
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
        
        graph_filename = f"Relio_maze_{JOB_ID}.png"
        graph_path = os.path.join(IMAGES_DIR, graph_filename)
        plt.savefig(graph_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"    ✔ Graph saved as '{graph_filename}'")
        return graph_filename
    
    except Exception as e:
        print(f"    ⚠️ Graph generation error: {e}")
        plt.close('all')
        return None

# =========================
# REPORT GENERATION
# =========================
def generate_full_report(compound, outcome, gene_evidence, pathway_map, paper_ids, mode_name):
    """Generate report as dictionary and text file"""
    print("\n[5] 📝 Generating Report...")
    fingerprint = compute_fingerprint(gene_evidence, pathway_map)
    
    conflicts = []
    if "FULL" in mode_name:
        conflicts = analyze_contradictions(gene_evidence)
    
    gene_to_paths = defaultdict(list)
    for p, genes in pathway_map.items():
        for g in genes:
            gene_to_paths[g].append(p)
    
    # Build report lines for text file
    lines = []
    lines.append("=" * 80)
    lines.append(f"RELIO V0.1 REPORT: {compound.upper()} + {outcome.upper()}")
    lines.append("=" * 80)
    lines.append(f"Mode             : {mode_name}")
    lines.append(f"Docs Analyzed    : {len(paper_ids)} (Strict Limit: 100)")
    lines.append(f"Net Direction    : {fingerprint['direction']}")
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
        if any(c['gene'] == g for c in conflicts):
            d_mode = "MIXED"
        
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
            if "PMC" not in pid_str and pid_str != "Unknown":
                pid_str = "PMC" + pid_str
            if "PMC" in pid_str:
                lines.append(f"- https://www.ncbi.nlm.nih.gov/pmc/articles/{pid_str}/")
        else:
            lines.append(f"- https://pubmed.ncbi.nlm.nih.gov/{pid_str}/")
    
    final_report_str = "\n".join(lines)
    
    # Save text report
    report_path = os.path.join(RESULTS_DIR, f"Relio_report_{JOB_ID}.txt")
    with open(report_path, "w") as f:
        f.write(final_report_str)
    print(f"    ✔ Report saved as 'Relio_report_{JOB_ID}.txt'")
    
    # Print to terminal
    print("\n" + "=" * 30 + " TERMINAL REPORT OUTPUT " + "=" * 30 + "\n")
    print(final_report_str)
    print("\n" + "=" * 80 + "\n")
    
    # Build JSON report for API
    gene_list = []
    for g, ev in top_genes_all:
        hits = len(ev)
        dirs = [e['direction'] for e in ev]
        d_mode = Counter(dirs).most_common(1)[0][0]
        if any(c['gene'] == g for c in conflicts):
            d_mode = "MIXED"
        
        my_paths = gene_to_paths.get(g, [])
        gene_list.append({
            "gene": g,
            "hits": hits,
            "direction": d_mode,
            "pathways": my_paths
        })
    
    source_links = []
    for pid in unique_ids:
        pid_str = str(pid)
        if mode_name.startswith("FULL"):
            if "PMC" not in pid_str and pid_str != "Unknown":
                pid_str = "PMC" + pid_str
            if "PMC" in pid_str:
                source_links.append({
                    "id": pid_str,
                    "url": f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pid_str}/"
                })
        else:
            source_links.append({
                "id": pid_str,
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pid_str}/"
            })
    
    json_report = {
        "job_id": JOB_ID,
        "compound": compound,
        "outcome": outcome,
        "mode": mode_name,
        "docs_analyzed": len(paper_ids),
        "fingerprint": fingerprint,
        "contradictions": conflicts if "FULL" in mode_name else [],
        "genes": gene_list,
        "sources": source_links,
        "images": {}
    }
    
    return json_report

# =========================
# MAIN EXECUTION
# =========================
def main():
    """Main execution function"""
    compound = os.environ.get("COMPOUND", "").strip()
    outcome = os.environ.get("OUTCOME", "").strip()
    mode = os.environ.get("MODE", "FAST").strip().upper()
    
    if not compound or not outcome:
        result = {
            "job_id": JOB_ID,
            "error": "Missing compound or outcome parameters",
            "status": "failed"
        }
        save_result(result)
        sys.exit(1)
    
    print("\n" + "=" * 80)
    print(f"🔬 RELIO V0.1 - GitHub Actions Runner")
    print(f"Compound: {compound}")
    print(f"Outcome:  {outcome}")
    print(f"Mode:     {mode}")
    print(f"Job ID:   {JOB_ID}")
    print("=" * 80)
    
    # Step 1: Build whitelist
    whitelist = build_gene_whitelist(outcome)
    
    if not whitelist:
        result = {
            "job_id": JOB_ID,
            "error": "Failed to build gene whitelist. Check outcome parameter.",
            "status": "failed"
        }
        save_result(result)
        sys.exit(1)
    
    # Step 2: Mine evidence
    if mode == "FAST":
        gene_evidence, paper_ids = mine_abstracts_fast(compound, outcome, whitelist, limit=100)
    else:
        gene_evidence, paper_ids = mine_pmc_full(compound, outcome, whitelist, limit=100)
    
    if not gene_evidence:
        result = {
            "job_id": JOB_ID,
            "compound": compound,
            "outcome": outcome,
            "mode": mode,
            "error": "No evidence found in database",
            "status": "no_results"
        }
        save_result(result)
        sys.exit(0)
    
    # Step 3: Build pathway map
    pathway_map = build_pathway_map(gene_evidence)
    
    # Step 4: Generate graph (if FULL mode)
    graph_filename = draw_graph(compound, gene_evidence, pathway_map, mode, outcome)
    
    # Step 5: Generate report
    json_report = generate_full_report(compound, outcome, gene_evidence, pathway_map, paper_ids, mode)
    
    # Add image reference if generated
    if graph_filename:
        json_report["images"]["maze"] = graph_filename
    
    json_report["status"] = "success"
    
    # Save JSON result
    save_result(json_report)
    
    print("\n✅ Analysis complete!")

def save_result(result):
    """Save result JSON to file"""
    result_path = os.path.join(RESULTS_DIR, f"{JOB_ID}.json")
    with open(result_path, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\n✔ Result saved: {result_path}")

if __name__ == "__main__":
    main()
