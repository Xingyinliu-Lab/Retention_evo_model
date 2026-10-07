#!/usr/bin/env python3
"""Exact paired AIC tests and fixed-Q ARD transition-route decomposition."""

import csv, io, itertools, json, math
from pathlib import Path

import numpy as np
from Bio import Phylo
from scipy.linalg import expm

ROOT=Path(__file__).resolve().parents[1]
FIT=ROOT/"results/MreB_ALP_composition_corrected_five_models_six_trees.tsv"
INPUT=ROOT/"input/six_tree_states.json"
OUT=ROOT/"results"
LEVELS=("M-only","M+ALP coexist","ALP-only")
INDEX={x:i for i,x in enumerate(LEVELS)}
GROUPS={"coexist_mediated":((0,1),(1,0),(1,2),(2,1)),"direct_endpoint":((0,2),(2,0))}

def read_tsv(path):
    with path.open(encoding="utf-8-sig",newline="") as h:return list(csv.DictReader(h,delimiter="\t"))
def write_tsv(path,rows):
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n");w.writeheader();w.writerows(rows)
def stationary(q):
    v,w=np.linalg.eig(q.T);x=np.abs(np.real(w[:,np.argmin(np.abs(v))]));return x/x.sum()
def parse_q(text):
    q=np.zeros((3,3))
    for token in text.split(";"):
        edge,value=token.rsplit(":",1);a,b=edge.split("->",1);q[INDEX[a],INDEX[b]]=float(value)
    np.fill_diagonal(q,-q.sum(axis=1));return q

class Engine:
    def __init__(self,tree,states,prior):
        self.nodes=list(tree.find_clades(order="postorder"));idx={id(n):i for i,n in enumerate(self.nodes)}
        self.children=[[idx[id(c)] for c in n.clades] for n in self.nodes]
        self.names=[n.name if n.is_terminal() else None for n in self.nodes]
        self.states=states;self.prior=prior;self.branch_count=sum(map(len,self.children))
    def logl(self,matrix):
        p=expm(matrix);likes=[np.ones(3) for _ in self.nodes];scale=0.0
        for i,n in enumerate(self.nodes):
            if n.is_terminal():
                v=np.zeros(3);v[self.states[self.names[i]]]=1;likes[i]=v;continue
            v=np.ones(3)
            for child in self.children[i]:v*=p@likes[child]
            z=v.sum()
            if z<=0 or not np.isfinite(z):raise FloatingPointError(z)
            likes[i]=v/z;scale+=math.log(z)
        root=float(self.prior@likes[-1]);return math.log(root)+scale
def tilted(q,edges,theta):
    m=q.copy();factor=math.exp(theta)
    for a,b in edges:m[a,b]=q[a,b]*factor
    np.fill_diagonal(m,np.diag(q));return m
def moments(engine,q,edges,h=5e-4):
    ll0=engine.logl(q);d={}
    for k in (-2,-1,1,2):d[k]=engine.logl(tilted(q,edges,k*h))-ll0
    mean=(d[-2]-8*d[-1]+8*d[1]-d[2])/(12*h)
    var=max(0.0,(-d[2]+16*d[1]+16*d[-1]-d[-2])/(12*h*h))
    return mean,var,ll0
def exact_signed_rank(values):
    values=[x for x in values if abs(x)>1e-12];n=len(values)
    order=sorted(range(n),key=lambda i:abs(values[i]));ranks=[0.0]*n;j=0
    while j<n:
        k=j+1
        while k<n and abs(abs(values[order[k]])-abs(values[order[j]]))<1e-12:k+=1
        rank=(j+1+k)/2
        for z in order[j:k]:ranks[z]=rank
        j=k
    obs=sum(r for r,x in zip(ranks,values) if x>0);total=sum(ranks)
    stats=[sum(r for r,s in zip(ranks,signs) if s) for signs in itertools.product((0,1),repeat=n)]
    return sum(abs(x-total/2)>=abs(obs-total/2)-1e-12 for x in stats)/len(stats)

def main():
    rows=read_tsv(FIT);trees=sorted({r["tree"] for r in rows});bench="adjacency_coexist_middle"
    lookup={(r["tree"],r["model"]):float(r["AIC"]) for r in rows};paired=[]
    for model in sorted({r["model"] for r in rows}-{bench}):
        diffs=[lookup[(t,model)]-lookup[(t,bench)] for t in trees]
        paired.append({"benchmark":bench,"comparison_model":model,"n_trees":len(trees),
                       "retention_better_trees":sum(x>0 for x in diffs),"comparison_better_trees":sum(x<0 for x in diffs),
                       "median_AIC_other_minus_retention":float(np.median(diffs)),"exact_p_two_sided":exact_signed_rank(diffs)})
    write_tsv(OUT/"paired_tests_vs_sequential_retention.tsv",paired)

    payload=json.loads(INPUT.read_text(encoding="utf-8"));items={x["name"]:x for x in payload["trees"]}
    summaries=[];groups=[]
    for fit in [r for r in rows if r["model"]=="ARD_unrestricted"]:
        item=items[fit["tree"]];tree=Phylo.read(io.StringIO(item["newick"]),"newick")
        states={k:INDEX[v] for k,v in item["states"].items()};q=parse_q(fit["rates"]);engine=Engine(tree,states,stationary(q))
        values={}
        for group,edges in GROUPS.items():
            mean,var,ll=moments(engine,q,edges);values[group]=(mean,var)
            groups.append({"tree":fit["tree"],"transition_group":group,"posterior_mean":mean,"posterior_variance":var,
                           "posterior_sd":math.sqrt(var),"log_likelihood_recomputed":ll,
                           "fit_log_likelihood":fit["log_likelihood"],"log_likelihood_abs_difference":abs(ll-float(fit["log_likelihood"]))})
        adjacent,va=values["coexist_mediated"];direct,vd=values["direct_endpoint"];total=adjacent+direct
        summaries.append({"tree":fit["tree"],"coexist_mediated_mean":adjacent,"direct_endpoint_mean":direct,
                          "coexist_mediated_fraction":adjacent/total,"direct_endpoint_fraction":direct/total,
                          "expected_total":total,"expected_total_per_branch":total/engine.branch_count})
    write_tsv(OUT/"MreB_ALP_ARD_transition_group_posterior_moments_six_trees.tsv",groups)
    write_tsv(OUT/"MreB_ALP_ARD_transition_contribution_summary_six_trees.tsv",summaries)
    if max(float(r["log_likelihood_abs_difference"]) for r in groups)>1e-4:raise SystemExit("ARD likelihood reproduction failed")
    print(json.dumps({"status":"complete","paired_tests":len(paired),"ARD_trees":len(summaries)}))
if __name__=="__main__":main()
