import sys, os, json, shutil, time
from python_runner.litmap import *

# Print exact Colab banner
print("""
LITMAP V0.1
===========
- TOOL NAME: LitMap V0.1
- FAST MODE: Max 100 Abstracts. Speed focused. (No Graph).
- FULL MODE: Max 100 Full Text. Maze Graph with LEGENDS. Contradiction Analysis.
- OUTPUT: Terminal Print + Text File + PNG Graph.
""")

compound, outcome, mode, job_id = sys.argv[1:]

BASE = os.path.dirname(__file__)
IMG_DIR = os.path.join(BASE, "images")
RES_DIR = os.path.join(BASE, "results")
os.makedirs(IMG_DIR, exist_ok=True)
os.makedirs(RES_DIR, exist_ok=True)

print("RELIO JOB:", compound, outcome, mode, job_id)
print()

whitelist = build_gene_whitelist(outcome)

mode_name = None
if whitelist:
    if mode == "FAST":
        ev, ids = mine_abstracts_fast(compound, outcome, whitelist, limit=100)
        mode_name = "FAST (Abstracts)"
    else:
        ev, ids = mine_pmc_full(compound, outcome, whitelist, limit=100)
        mode_name = "FULL (PMC Full Text)"

    if ev:
        pmap = build_pathway_map(ev)
        draw_graph(compound, ev, pmap, mode_name, outcome)
        generate_full_report(compound, outcome, ev, pmap, ids, mode_name)
    else:
        print("❌ No evidence found.")
else:
    print("❌ Whitelist failed.")

# Move outputs with job_id (exact same logic)
src_img = os.path.join(IMG_DIR, "litmap_maze.png")
dst_img = os.path.join(IMG_DIR, f"{job_id}_graph.png")
if os.path.exists(src_img):
    shutil.move(src_img, dst_img)

src_rep = os.path.join(RES_DIR, "litmap_report.txt")
dst_rep = os.path.join(RES_DIR, f"{job_id}_report.txt")
if os.path.exists(src_rep):
    shutil.move(src_rep, dst_rep)

# Generate JSON result (unchanged)
result = {
    "compound": compound,
    "outcome": outcome,
    "mode": mode,
    "mode_name": mode_name,
    "genes": dict(ev) if 'ev' in locals() else {},
    "pathways": {k: list(v) for k, v in pmap.items()} if 'pmap' in locals() else {},
    "documents": ids if 'ids' in locals() else [],
    "images": {"graph": f"{job_id}_graph.png" if os.path.exists(dst_img) else None},
    "job_id": job_id
}

with open(os.path.join(RES_DIR, f"{job_id}.json"), "w") as f:
    json.dump(result, f, indent=2)

print("RELIO JOB COMPLETE:", job_id)
