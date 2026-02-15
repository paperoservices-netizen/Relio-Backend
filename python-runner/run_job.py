import sys,json,os
from relio import run_relio

compound,outcome,mode,job_id=sys.argv[1:5]
res=run_relio(compound,outcome,mode,job_id)

os.makedirs("python-runner/results",exist_ok=True)
with open(f"python-runner/results/{job_id}.json","w") as f:
    json.dump(res,f,indent=2)
