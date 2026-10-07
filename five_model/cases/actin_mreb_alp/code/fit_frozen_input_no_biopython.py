#!/usr/bin/env python3
"""Fit the frozen six-tree MreB/ALP input; requires only NumPy and SciPy."""

import csv, json, math, multiprocessing as mp, os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
LEVELS = ("M-only", "M+ALP coexist", "ALP-only")
CODE = {state: index for index, state in enumerate(LEVELS)}
MODELS = [
    ("ARD_unrestricted", "stationary", tuple((a,b) for a in range(3) for b in range(3) if a != b)),
    ("adjacency_coexist_middle", "stationary", ((0,1),(1,0),(1,2),(2,1))),
    ("adjacency_ALP_middle", "stationary", ((0,2),(2,0),(2,1),(1,2))),
    ("adjacency_M_middle", "stationary", ((2,0),(0,2),(0,1),(1,0))),
    ("ancestral_coexistence_strict", 1, ((1,0),(1,2))),
]

class Node:
    def __init__(self):
        self.name, self.children = "", []
    def terminal(self):
        return not self.children

def parse_newick(text):
    text=text.strip().rstrip(";"); pos=0
    def parse():
        nonlocal pos
        node=Node()
        if text[pos:pos+1] == "(":
            pos += 1
            while True:
                node.children.append(parse())
                if text[pos:pos+1] == ",": pos += 1; continue
                if text[pos:pos+1] == ")": pos += 1; break
        start=pos
        while pos < len(text) and text[pos] not in ",()": pos += 1
        token=text[start:pos].strip()
        node.name=token.rsplit(":",1)[0] if ":" in token else token
        return node
    return parse()

def stationary(q):
    values,vectors=np.linalg.eig(q.T)
    vector=np.abs(np.real(vectors[:,np.argmin(np.abs(values))]))
    return vector/vector.sum()

def fit(task):
    item,model,root_mode,edges=task
    root=parse_newick(item["newick"]); tip_states={k:CODE[v] for k,v in item["states"].items()}
    nodes=[]
    def walk(node):
        for child in node.children: walk(child)
        nodes.append(node)
    walk(root); index={id(node):i for i,node in enumerate(nodes)}
    children=[[index[id(child)] for child in node.children] for node in nodes]
    names=[node.name if node.terminal() else None for node in nodes]
    def objective(log_rates):
        q=np.zeros((3,3))
        for value,(a,b) in zip(np.exp(log_rates),edges): q[a,b]=value
        np.fill_diagonal(q,-q.sum(axis=1)); transition=expm(q)
        likes=[np.ones(3) for _ in nodes]; logscale=0.0
        for i,node in enumerate(nodes):
            if node.terminal():
                vector=np.zeros(3); vector[tip_states[names[i]]]=1; likes[i]=vector; continue
            vector=np.ones(3)
            for child_i in children[i]: vector *= transition @ likes[child_i]
            total=vector.sum()
            if total <= 0 or not np.isfinite(total): return 1e100
            likes[i]=vector/total; logscale += math.log(total)
        prior=stationary(q) if root_mode == "stationary" else np.eye(3)[root_mode]
        value=float(prior @ likes[-1])
        return -(math.log(value)+logscale) if value > 0 and np.isfinite(value) else 1e100
    k=len(edges)
    starts=[np.zeros(k),np.full(k,-1.0),np.full(k,1.0),np.linspace(-2,2,k),np.linspace(2,-2,k)]
    fits=[minimize(objective,start,method="L-BFGS-B",bounds=[(-10,10)]*k,
                   options={"maxiter":400,"ftol":1e-11,"gtol":1e-7}) for start in starts]
    converged=[result for result in fits if result.success and np.isfinite(result.fun)]
    best=min(converged or fits,key=lambda result:result.fun)
    counts={state:sum(value==state for value in item["states"].values()) for state in LEVELS}
    return {"tree":item["name"],"model":model,
            "root_mode":LEVELS[root_mode] if isinstance(root_mode,int) else root_mode,
            "branch_scale":"unit_topology","k":k,"log_likelihood":-float(best.fun),
            "AIC":2*k+2*float(best.fun),"success":str(bool(best.success)).upper(),
            "gradient_max":float(np.max(np.abs(best.jac))) if best.jac is not None else "",
            "rates":";".join(f"{LEVELS[a]}->{LEVELS[b]}:{math.exp(x):.10g}" for x,(a,b) in zip(best.x,edges)),
            "mapped_MAGs":len(tip_states),"root_projection_mismatch":item["root_projection_mismatch"],
            "n_M_only":counts["M-only"],"n_coexist":counts["M+ALP coexist"],"n_ALP_only":counts["ALP-only"]}

def main():
    payload=json.loads((ROOT/"input/six_tree_states.json").read_text(encoding="utf-8"))
    tasks=[(item,*model) for item in payload["trees"] for model in MODELS]
    with mp.Pool(processes=len(tasks)) as pool: rows=pool.map(fit,tasks)
    for tree in {row["tree"] for row in rows}:
        block=[row for row in rows if row["tree"]==tree]; minimum=min(row["AIC"] for row in block)
        for row in block: row["delta_AIC"]=row["AIC"]-minimum
    out=ROOT/"results"; out.mkdir(parents=True,exist_ok=True)
    with (out/"MreB_ALP_composition_corrected_five_models_six_trees.tsv").open("w",encoding="utf-8",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    (ROOT/"FIT_COMPLETE.marker").write_text("status=complete\nmodels=5\ntrees=6\nbranch_scale=unit_topology\nrunner=no_biopython\n",encoding="utf-8")
    print(json.dumps({"status":"complete","rows":len(rows),"processes":len(tasks)}))

if __name__ == "__main__": main()
