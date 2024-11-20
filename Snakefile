import json
import itertools

configfile: "config.yaml"
R = config["R"]

# WILDCARDS --------------------------------------------------------------------
# tsb = simulated type/state/batch effect
# sim = simulated/synthetic dataset identifier
# dat = experimental/real dataset identifier
# sco = gene-level type-/state-specificity score
# sel = feature selection method for reprocessing
# sta = downstream evaluation statistic
# das = differential abundance/state method
# ------------------------------------------------------------------------------

# magnitude of type, state & batch effect
T = list(range(0, 120, 20)); S = list(range(0, 120, 20)); B = [0]
SIM = ["t{},s{},b{}".format(t,s,b) for t in T for s in S for b in B]
NUM = list(range(1, 10, 1))

DAT = glob_wildcards("code/00-get_dat-{x}.R").x
SCO = glob_wildcards("code/02-sco-{x}.R").x
SEL = glob_wildcards("code/03-sel_sim-{x}.R").x
STA = glob_wildcards("code/05-sta-{x}.R").x
DAS = glob_wildcards("code/06-das-{x}.R").x

# specify scores required for selections,
# defaulting to 'random' if left unspecified
sco_by_sel = {
    "HVG": ["HVG"],
    "Fstat": ["tF"],
    "tPVE_sPVE": ["PVE"],
    "tF_sPVE": ["tF", "PVE"],
    "tF_sPBDS": ["tF", "sPBDS"]}
for sel in SEL:
    if sel not in sco_by_sel.keys():
        sco_by_sel[sel] = ["random"]

def sco_by_sim(wildcards):
    return expand(sim_out+"sco-{sim},{sco}.rds", 
        sim=wildcards.sim, sco=sco_by_sel[wildcards.sel])

def sco_by_dat(wildcards):
    return expand(dat_out+"sco-{dat},{sco}.rds",
        dat=wildcards.dat, sco=sco_by_sel[wildcards.sel])


# output directories
sim_dat = "data/sim/"; dat_dat = "data/dat/"
sim_out = "outs/sim/"; dat_out = "outs/dat/"

# simulation
sim = expand(sim_dat+"00-raw/t{t},s{s},b{b}.rds", t=T, s=S, b=B)
sim_pro = expand(sim_dat+"01-pro/{sim}.rds", sim=SIM)
sim_pro_rd = expand(sim_dat+"01-pro/{sim}-rd.rds", sim=SIM)
sim_pro_cd = expand(sim_dat+"01-pro/{sim}-cd.rds", sim=SIM)
sim_sco = expand(sim_out+"sco-{sim},{sco}.rds", sim=SIM, sco=SCO)
sim_sel = expand(sim_out+"sel-{sim},{sel}.rds", sim=SIM, sel=SEL)
sim_rep = expand(sim_dat+"02-rep/{sim},{sel}.rds", sim=SIM, sel=SEL)
sim_rep_cd = expand(sim_dat+"02-rep/{sim},{sel}-cd.rds", sim=SIM, sel=SEL)
sim_sta = expand(sim_out+"sta-{sim},{sel},{sta}.rds", sim=SIM, sel=SEL, sta=STA)
sim_das = expand(sim_out+"das-{sim},{sel},{das}.rds", sim=SIM, sel=SEL, das=DAS)

# for real data...
# exclude simulation parameter-based selection
SEL = glob_wildcards("code/03-sel_dat-{x}.R").x
# exclude ground-truth based evaluation statistics
STA = [x for x in STA if "F1" not in x]

# application
dat = expand(dat_dat+"00-raw/{dat}.rds", dat=DAT)
dat_pro = expand(dat_dat+"01-pro/{dat}.rds", dat=DAT)
dat_pro_rd = expand(dat_dat+"01-pro/{dat}-rd.rds", dat=DAT)
dat_pro_cd = expand(dat_dat+"01-pro/{dat}-cd.rds", dat=DAT)
dat_sco = expand(dat_out+"sco-{dat},{sco}.rds", dat=DAT, sco=SCO)
dat_sel = expand(dat_out+"sel-{dat},{sel}.rds", dat=DAT, sel=SEL)
dat_rep = expand(dat_dat+"02-rep/{dat},{sel}.rds", dat=DAT, sel=SEL)
dat_rep_cd = expand(dat_dat+"02-rep/{dat},{sel}-cd.rds", dat=DAT, sel=SEL)
dat_sta = expand(dat_out+"sta-{dat},{sel},{sta}.rds", dat=DAT, sel=SEL, sta=STA)
dat_das = expand(dat_out+"das-{dat},{sel},{das}.rds", dat=DAT, sel=SEL, das=DAS)
dat_eva = expand(dat_out+"eva-{dat},{sel},{sta},{num}.rds", dat=DAT, sel=SEL, sta=STA, num=NUM)

sim_res = {
    "sim": sim,
    "pro": sim_pro,
    "rd1": sim_pro_rd,
    "cd1": sim_pro_cd,
    "sco": sim_sco,
    "sel": sim_sel,
    "rep": sim_rep,
    "cd2": sim_rep_cd,
    "sta": sim_sta,
    "das": sim_das}

dat_res = {
    "dat": dat,
    "pro": dat_pro,
    "sco": dat_sco,
    "sel": dat_sel,
    "rep": dat_rep,
    "cd2": dat_rep_cd,
    "sta": dat_sta,
    "eva": dat_eva,
    "das": dat_das}

# visualization
VAL = sim_res.keys()
WAL = dat_res.keys()

plt = []
for val in VAL:
    x = glob_wildcards("code/08-plt_"+val+"-{x}.R").x
    plt += expand("plts/sim/{val}-{plt}.pdf", val=val, plt=x)

qlt = []
for wal in WAL:
    x = glob_wildcards("code/08-qlt_"+wal+"-{x}.R").x
    qlt += expand("plts/dat/{wal}-{qlt}.pdf", wal=wal, qlt=x)

# SETUP ========================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
        "session_info.txt",
        [x for x in sim_res.values()], plt,
        [x for x in dat_res.values()], qlt

rule session_info:
    priority: 100
    input:  "code/10-session_info.R"
    output: "session_info.txt"
    log:    "logs/session_info.Rout"
    shell:  '''
    {R} CMD BATCH --no-restore --no-save\
    "--args {output}" {input} {log}'''

# SIMULATION ===================================================================

rule get_sim:
    priority: 99
    input:  "code/00-get_sim.R",
    output: sim_dat+"00-raw/t{t},s{s},b{b}.rds"
    log:    "logs/get_sim-t{t},s{s},b{b}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} res={output}" {input[0]} {log}'''

# processing
rule pro_sim:
    priority: 98
    input:  "code/01-pro_sim.R",
            sim_dat+"00-raw/{sim}.rds"
    output: sim_dat+"01-pro/{sim}.rds",
            sim_dat+"01-pro/{sim}-rd.rds",
            sim_dat+"01-pro/{sim}-cd.rds"
    log:    "logs/pro_sim-{sim}.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {output[0]} {output[1]} {output[2]}" {input[0]} {log}'''

# scoring
rule sco_sim:
    priority: 97
    input:  "code/02-sco.R",
            "code/02-sco-{sco}.R",
            rules.pro_sim.output[0]
    output: sim_out+"sco-{sim},{sco}.rds"
    log:    "logs/sco_sim-{sim},{sco}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

def sco_by_sim(wildcards):
    return expand(
        sim_out+"sco-{sim},{sco}.rds", 
        sim=wildcards.sim, sco=SCO)

# selection
rule sel_sim:
    priority: 96
    input:  "code/03-sel.R",
            "code/03-sel_sim-{sel}.R",
            x = sco_by_sim
    params: lambda wc, input: ";".join(input.x)
    output: sim_out+"sel-{sim},{sel}.rds"
    log:    "logs/sel_sim-{sim},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {params} {output[0]}" {input[0]} {log}'''

# reprocessing
rule rep_sim:
    priority: 95
    input:  "code/04-rep.R",
            rules.pro_sim.output[0],
            rules.sel_sim.output
    output: sim_dat+"02-rep/{sim},{sel}.rds",
            sim_dat+"02-rep/{sim},{sel}-cd.rds"
    log:    "logs/rep_sim-{sim},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output[0]} {output[1]}" {input[0]} {log}'''

# evaluation
rule sta_sim:
    priority: 94
    input:  "code/05-sta.R",
            "code/05-sta-{sta}.R",
            rules.rep_sim.output
    output: sim_out+"sta-{sim},{sel},{sta}.rds"
    log:    "logs/sta_sim-{sim},{sel},{sta}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# differential
rule das_sim:
    priority: 93
    input:  "code/06-das.R",
            "code/06-das-{das}.R",
            rules.rep_sim.output
    output: sim_out+"das-{sim},{sel},{das}.rds"
    log:    "logs/das_sim-{sim},{sel},{das}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# APPLICATION ==================================================================

rule get_dat:
    priority: 49
    input:  "code/00-get_dat.R",
            "code/00-get_dat-{dat}.R"
    output: dat_dat+"00-raw/{dat}.rds"
    log:	"logs/get_dat-{dat}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {input[1]} {output}" {input[0]} {log}'''    

# processing
rule pro_dat:
    priority: 48
    input:  "code/01-pro_dat.R",
            rules.get_dat.output
    output: dat_dat+"01-pro/{dat}.rds",
            dat_dat+"01-pro/{dat}-rd.rds",
            dat_dat+"01-pro/{dat}-cd.rds"
    log:    "logs/pro_dat-{dat}.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {output[0]} {output[1]} {output[2]}" {input[0]} {log}'''

# scoring
rule sco_dat:
    priority: 47
    input:  "code/02-sco.R",
            "code/02-sco-{sco}.R",
            rules.pro_dat.output[0]
    output: dat_out+"sco-{dat},{sco}.rds"
    log:    "logs/sco_dat-{dat},{sco}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

def sco_by_dat(wildcards):
    return expand(
        dat_out+"sco-{dat},{sco}.rds",
        dat=wildcards.dat, sco=SCO)

# selection
rule sel_dat:
    priority: 46
    input:  "code/03-sel.R",
            "code/03-sel_dat-{sel}.R",
            x = sco_by_dat
    params: lambda wc, input: ";".join(input.x)
    output: dat_out+"sel-{dat},{sel}.rds"
    log:    "logs/sel_dat-{dat},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {params} {output[0]}" {input[0]} {log}'''


# reprocessing
rule rep_dat:
    priority: 45
    input:  "code/04-rep.R",
            rules.pro_dat.output[0],
            rules.sel_dat.output
    output: dat_dat+"02-rep/{dat},{sel}.rds",
            dat_dat+"02-rep/{dat},{sel}-cd.rds"
    log:    "logs/rep_dat-{dat},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output[0]} {output[1]}" {input[0]} {log}'''

# determine the optimal number of genes
rule eva_dat:
    priority: 44
    input:  "code/07-eva.R",
            "code/05-sta-{sta}.R",
            rules.pro_dat.output[0],
            x = sco_by_dat
    params: lambda wc, input: ";".join(input.x)
    output: dat_out+"eva-{dat},{sel},{sta},{num}.rds"
    log:    "logs/eva_dat-{dat},{sel},{sta},{num}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {params} {output[0]}" {input[0]} {log}'''


# evaluation
rule sta_dat:
    priority: 43
    input:  "code/05-sta.R",
            "code/05-sta-{sta}.R",
            rules.rep_dat.output
    output: dat_out+"sta-{dat},{sel},{sta}.rds"
    log:    "logs/sta_dat-{dat},{sel},{sta}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''


# differential
rule das_dat:
   priority: 42
   input:  "code/06-das.R",
           "code/06-das-{das}.R",
           rules.rep_dat.output
   output: dat_out+"das-{dat},{sel},{das}.rds"
   log:    "logs/das_dat-{dat},{sel},{das}.Rout"
   shell: '''
       {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
       {input[1]} {input[2]} {output}" {input[0]} {log}'''
    

# VISUALIZATION ========================================================

for val in VAL:
    rule:
        priority: 90
        input:  expand("code/08-plt_{val}-{{plt}}.R", val=val), x=sim_res[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/sim/{val}-{{plt}}.pdf", val=val)
        log:    expand("logs/plt_{val}-{{plt}}.Rout", val=val)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''

for wal in WAL:
    rule:
        priority: 40
        input:  expand("code/08-qlt_{wal}-{{qlt}}.R", wal=wal), x=dat_res[wal]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/dat/{wal}-{{qlt}}.pdf", wal=wal)
        log:    expand("logs/qlt_{wal}-{{qlt}}.Rout", wal=wal)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''
