import sys, json, os
from litmap_engine import run_litmap_job

input_path = sys.argv[1]
job_id = sys.argv[2]

with open(input_path) as f:
    data = json.load(f)

# Create job folders
result_dir = f"python-runner/results/{job_id}"
image_dir = f"python-runner/images/{job_id}"
os.makedirs(result_dir, exist_ok=True)
os.makedirs(image_dir, exist_ok=True)

run_litmap_job(data, job_id, result_dir, image_dir)
