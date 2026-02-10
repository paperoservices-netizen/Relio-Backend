import sys, json, os
from relio_engine import run_litmap

compound = sys.argv[1]
outcome  = sys.argv[2]
mode     = sys.argv[3].upper()
job_id   = sys.argv[4]

print("🧬 Relio job:", job_id)

os.makedirs("python-runner/results/images", exist_ok=True)

result = run_litmap(
    compound=compound,
    outcome=outcome,
    mode=mode,
    job_id=job_id
)

with open(f"python-runner/results/{job_id}.json", "w") as f:
    json.dump(result, f, indent=2)