import sys, os, json, shutil
from litmap import *

compound, outcome, mode, job_id = sys.argv[1:]

os.makedirs("python_runner/images",exist_ok=True)
os.makedirs("python_runner/results",exist_ok=True)

whitelist = build_gene_whitelist(outcome)

if mode=="FAST":
    genes,docs = mine_abstracts_fast(compound,outcome,whitelist)
else:
    genes,docs = mine_pmc_full(compound,outcome,whitelist)

pathways = build_pathway_map(genes)
draw_graph(compound,genes,pathways,mode,outcome)
generate_full_report(compound,outcome,genes,pathways,docs,mode)

if os.path.exists("python_runner/images/litmap_maze.png"):
    shutil.move("python_runner/images/litmap_maze.png",f"python_runner/images/{job_id}_graph.png")

if os.path.exists("python_runner/results/litmap_report.txt"):
    shutil.move("python_runner/results/litmap_report.txt",f"python_runner/results/{job_id}_report.txt")

result={
 "compound":compound,
 "outcome":outcome,
 "mode":mode,
 "genes":genes,
 "pathways":{k:list(v) for k,v in pathways.items()},
 "documents":docs,
 "fingerprint":compute_fingerprint(genes,pathways),
 "contradictions":analyze_contradictions(genes),
 "images":{"graph":f"{job_id}_graph.png"}
}

with open(f"python_runner/results/{job_id}.json","w") as f:
    json.dump(result,f,indent=2)

print("RELIO JOB DONE",job_id)
