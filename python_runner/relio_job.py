import sys, os, json, shutil
from python_runner.litmap import *

compound, outcome, mode, job_id = sys.argv[1:]

BASE = os.path.dirname(__file__)
IMG_DIR = os.path.join(BASE, "images")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(IMG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

print("RELIO JOB:", compound, outcome, mode, job_id)

whitelist = build_gene_whitelist(outcome)

if mode == "FAST":
    genes, docs = mine_abstracts_fast(compound, outcome, whitelist)
else:
    genes, docs = mine_pmc_full(compound, outcome, whitelist)

pathways = build_pathway_map(genes)

draw_graph(compound, genes, pathways, mode, outcome)
generate_full_report(compound, outcome, genes, pathways, docs, mode)

# Move outputs
src_img = os.path.join(IMG_DIR, "litmap_maze.png")
dst_img = os.path.join(IMG_DIR, f"{job_id}_graph.png")
if os.path.exists(src_img):
    shutil.move(src_img, dst_img)

src_rep = os.path.join(RES_DIR, "litmap_report.txt")
dst_rep = os.path.join(RES_DIR, f"{job_id}_report.txt")
if os.path.exists(src_rep):
    shutil.move(src_rep, dst_rep)

result = {
    "compound": compound,
    "outcome": outcome,
    "mode": mode,
    "genes": genes,
    "pathways": {k: list(v) for k, v in pathways.items()},
    "documents": docs,
    "fingerprint": compute_fingerprint(genes, pathways),
    "contradictions": analyze_contradictions(genes),
    "images": {
        "graph": f"{job_id}_graph.png"
    }
}

with open(os.path.join(RES_DIR, f"{job_id}.json"), "w") as f:
    json.dump(result, f, indent=2)

print("RELIO JOB COMPLETE:", job_id)
