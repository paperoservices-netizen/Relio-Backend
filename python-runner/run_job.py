import sys
from relio_engine import run_relio_job

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python run_job.py <compound> <outcome> <mode> <job_id>")
        sys.exit(1)

    compound = sys.argv[1]
    outcome = sys.argv[2]
    mode = sys.argv[3]
    job_id = sys.argv[4]

    print(f"\n🔬 Running job {job_id} with {compound} + {outcome} ({mode})")
    run_relio_job(compound, outcome, mode, job_id)
