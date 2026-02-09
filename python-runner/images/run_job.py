import sys, os
from litmap_engine import run_litmap_job

compound = sys.argv[1]
outcome = sys.argv[2]
mode = sys.argv[3].upper()
job_id = sys.argv[4]  # unique per run

result_dir = f"python-runner/results/{job_id}"
image_dir = f"python-runner/images/{job_id}"
os.makedirs(result_dir, exist_ok=True)
os.makedirs(image_dir, exist_ok=True)

data = {"compound": compound, "outcome": outcome, "mode": mode}

run_litmap_job(data, job_id, result_dir, image_dir)
