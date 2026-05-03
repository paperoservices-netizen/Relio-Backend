import sys
import os
import json
import shutil
from relio import *  # ensure relio.py is in the same directory or in your PYTHONPATH

print("""
RELIO V0.1
===========
- TOOL NAME: Relio V0.1
- FAST MODE: Max 100 Abstracts. Speed focused. (No Graph).
- FULL MODE: Max 100 Full Text. Maze Graph with LEGENDS. Contradiction Analysis.
- OUTPUT: Terminal Print + Text File + PNG Graph.
""")

# Validate arguments
if len(sys.argv) != 5:
    print("❌ ERROR: Expected 4 arguments: compound outcome mode job_id")
    print(f"Received: {len(sys.argv)-1} arguments")
    sys.exit(1)

compound, outcome, mode, job_id = sys.argv[1:]
print("RELIO JOB:", compound, outcome, mode, job_id)
print()

# Define directories
BASE = os.path.dirname(__file__)
IMG_DIR = os.path.join(BASE, "images")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(IMG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

# Run the main analysis as in your script
# Build whitelist
whitelist = build_gene_whitelist(outcome)
if not whitelist:
    print("❌ Whitelist build failed.")
    sys.exit(1)

# Run mining based on mode
if mode.upper() == "FAST":
    ev, ids = mine_abstracts_fast(compound, outcome, whitelist, limit=100)
    mode_name = "FAST (Abstracts)"
else:
    ev, ids = mine_pmc_full(compound, outcome, whitelist, limit=100)
    mode_name = "FULL (PMC Full Text)"

if not ev:
    print("❌ No evidence found.")
    sys.exit(0)

# Build pathway map
pmap = build_pathway_map(ev)

# Generate graph image
draw_graph(compound, ev, pmap, mode_name, outcome)

# Generate report
generate_full_report(compound, outcome, ev, pmap, ids, mode_name)

# Move generated files to the designated locations
src_img = os.path.join(IMG_DIR, "litmap_maze.png")
dst_img = os.path.join(IMG_DIR, f"{job_id}_graph.png")
if os.path.exists(src_img):
    shutil.move(src_img, dst_img)

src_rep = os.path.join(IMG_DIR, "litmap_report.txt")
dst_rep = os.path.join(RES_DIR, f"{job_id}_report.txt")
if os.path.exists(src_rep):
    shutil.move(src_rep, dst_rep)

# Prepare JSON output
pathways_serializable = {p: list(genes) for p, genes in pmap.items()}
result = {
    "compound": compound,
    "outcome": outcome,
    "mode": mode,
    "mode_name": mode_name,
    "genes": {gene: [{"id": h["id"], "direction": h["direction"]} for h in hits] for gene, hits in ev.items()},
    "pathways": pathways_serializable,
    "documents": list(set(ids)),
    "images": {"graph": f"{job_id}_graph.png" if os.path.exists(dst_img) else None},
    "job_id": job_id,
    "success": bool(ev),
    "total_genes": len(ev),
    "total_documents": len(ids),
    "total_pathways": len(pmap),
}

json_path = os.path.join(RES_DIR, f"{job_id}.json")
with open(json_path, "w") as f:
    json.dump(result, f, indent=2)

print("RELIO JOB COMPLETE:", job_id)
print(f"JSON saved: {json_path}")
