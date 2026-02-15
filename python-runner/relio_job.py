import sys, json, os, shutil
from litmap import (
    build_gene_whitelist,
    mine_abstracts_fast,
    mine_pmc_full,
    build_pathway_map,
    draw_graph,
    generate_full_report,
    analyze_contradictions,
    compute_fingerprint
)

compound = sys.argv[1]
outcome  = sys.argv[2]
mode     = sys.argv[3]
job_id   = sys.argv[4]

os.makedirs("python-runner/results", exist_ok=True)
os.makedirs("python-runner/images", exist_ok=True)

# Run LitMap
whitelist = build_gene_whitelist(outcome)

if "FAST" in mode:
    gene_evidence, docs = mine_abstracts_fast(compound, outcome, whitelist)
    mode_name = "FAST"
else:
    gene_evidence, docs = mine_pmc_full(compound, outcome, whitelist)
    mode_name = "FULL"

pathways = build_pathway_map(gene_evidence)

# Generate graph + report
draw_graph(compound, gene_evidence, pathways, mode_name, outcome)
generate_full_report(compound, outcome, gene_evidence, pathways, docs, mode_name)

# Move graph
graph_name = f"{job_id}_graph.png"
if os.path.exists("litmap_maze.png"):
    shutil.move("litmap_maze.png", f"python-runner/images/{graph_name}")

# Move report
if os.path.exists("litmap_report.txt"):
    shutil.move("litmap_report.txt", f"python-runner/results/{job_id}_report.txt")

# Build JSON
fingerprint = compute_fingerprint(gene_evidence, pathways)
contradictions = analyze_contradictions(gene_evidence) if "FULL" in mode else []

result = {
  "compound": compound,
  "outcome": outcome,
  "mode": mode,
  "fingerprint": fingerprint,
  "genes": gene_evidence,
  "pathways": {k:list(v) for k,v in pathways.items()},
  "contradictions": contradictions,
  "documents": docs,
  "images": {
      "graph": graph_name if os.path.exists(f"python-runner/images/{graph_name}") else None
  }
}

with open(f"python-runner/results/{job_id}.json","w") as f:
    json.dump(result,f,indent=2)

print("✅ RELIO JOB COMPLETE:", job_id)
