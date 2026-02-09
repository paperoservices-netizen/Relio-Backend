import sys, os
from relio_engine import run_relio_job

# Command-line arguments: compound, outcome, mode, job_id
compound = sys.argv[1]
outcome = sys.argv[2]
mode = sys.argv[3].upper()
job_id = sys.argv[4]

# Create job-specific result directories
result_dir = f"python-runner/results/{job_id}"
image_dir = f"python-runner/images/{job_id}"
os.makedirs(result_dir, exist_ok=True)
os.makedirs(image_dir, exist_ok=True)

data = {"compound": compound, "outcome": outcome, "mode": mode}

run_relio_job(data, job_id, result_dir, image_dir)
