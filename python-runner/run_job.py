import sys
import os
from relio_engine import run_relio_job

if len(sys.argv) != 5:
    print("Usage: python run_job.py <compound> <outcome> <mode> <job_id>")
    sys.exit(1)

compound = sys.argv[1]
outcome = sys.argv[2]
mode = sys.argv[3]
job_id = sys.argv[4]

# Create job-specific folders
os.makedirs(f"python-runner/results", exist_ok=True)
os.makedirs(f"python-runner/images/{job_id}", exist_ok=True)

# Run the analysis
result = run_relio_job(compound, outcome, mode, job_id)

# Save JSON result
import json
with open(f"python-runner/results/{job_id}.json", "w") as f:
    json.dump(result, f, indent=2)

print(f"✔ Job {job_id} completed. JSON saved in results/, graphs in images/{job_id}/")
