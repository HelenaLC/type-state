import json
import itertools

configfile: "config.yaml"
R = config["R"]

SCO = glob_wildcards("code/02-sco-{x}.R").x
SEL = glob_wildcards("code/03-sel-{x}.R").x
STA = glob_wildcards("code/05-sta-{x}.R").x
DAS = glob_wildcards("code/06-das-{x}.R").x
EVA = glob_wildcards("code/07-eva-{x}.R").x



# magnitude of type, state & batch effect
T = list(range(0, 120, 20))
S = list(range(0, 120, 20))
B = [0]

SIM = ["t{},s{},b{}".format(t,s,b) for t in T for s in S for b in B]
VAL = ["cd", "sco", "sta", "eva", "das"]
PLT = {val:[plt for plt in glob_wildcards("code/08-plot_" + val + "-{x}.R").x] for val in VAL}

res_sim = expand(
    "data/00-sim/t{t},s{s},b{b}.rds", 
    t = T, s = S, b = B)

res_fil = expand("data/01-fil/{sim}.rds",    sim = SIM)
res_rd  = expand("data/01-fil/{sim}-rd.rds", sim = SIM)
res_cd  = expand("data/01-fil/{sim}-cd.rds", sim = SIM)

res_sco = expand(
    "outs/sco-{sim},{sco}.rds",
    sim = SIM, sco = SCO)
res_sel = expand(
    "outs/sel-{sim},{sel}.rds",
    sim = SIM, sel = SEL)
res_rep = expand(
    "data/02-rep/{sim},{sel}.rds", 
    sim = SIM, sel = SEL)
res_sta = expand(
    "outs/sta-{sim},{sel},{sta}.rds",
    sim = SIM, sel = SEL, sta = STA)
res_das = expand(
    "outs/das-{sim},{sel},{das}.rds",
    sim = SIM, sel = SEL, das = DAS)
res_eva = expand(
    "outs/eva-{sim},{sco},{eva}.rds",
    sim = SIM, sco = SCO, eva = EVA)

res_plt = list()
for val in PLT.keys():
    res_plt += expand("plts/{val}-{plt}.pdf", val = val, plt = PLT[val])

res = {
    "sim": res_sim,
    "fil": res_fil,
    "rd":  res_rd, 
    "cd":  res_cd,
    "sco": res_sco,
    "sel": res_sel,
    "rep": res_rep,
    "sta": res_sta,
    "plt": res_plt,
    "eva": res_eva,
    "das": res_das,
    "plt": res_plt}



# COLLECTION ===================================================================

def res_sco_by_sim(wildcards):
    return expand("outs/sco-{sim},{sco}.rds",
        sim = wildcards.sim, sco = SCO)

def res_sta_by_sco(wildcards):
    return expand("outs/sta-{sim},{sco},{sel},{sta}.rds",
        sim = SIM, sco = wildcards.sco, sel = SEL, sta = STA)

# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
        "session_info.txt",
        # simulation, pre- & re-processing
        res_sim, res_fil, res_rep, 
        # scoring & evaluation
        res_sco, res_sel, res_sta, #res_roc,
        # downstream
        res_das,
        # visualization
        res_plt

rule session_info:
    priority: 100
    input:  "code/10-session_info.R"
    output: "session_info.txt"
    log:    "logs/session_info.Rout"
    shell:  '''
    {R} CMD BATCH --no-restore --no-save\
    "--args {output}" {input} {log}'''

# PREPROCESSING ================================================================

# synthetic data generation
rule sim_data:
    priority: 99
    input:  "code/00-sim_data.R",
    output: "data/00-sim/t{t},s{s},b{b}.rds"
    log:	"logs/sim_data-t{t},s{s},b{b}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        wcs={wildcards} res={output}" {input[0]} {log}'''

# filtering & preprocessing
rule fil_data:
    priority: 98
    input:  "code/01-fil_data.R",
            "data/00-sim/{sim}.rds"
    output: "data/01-fil/{sim}.rds",
            "data/01-fil/{sim}-rd.rds",
            "data/01-fil/{sim}-cd.rds"
    log:    "logs/fil_data-{sim}.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        sim={input[1]} fil={output[0]} rd={output[1]} cd={output[2]}" {input[0]} {log}'''

# ANALYSIS =====================================================================

# calculate variability scores
rule calc_sco:
    priority: 97
    input:  "code/02-sco.R",
            "code/02-sco-{sco}.R",
            "data/01-fil/{sim}.rds"
    output: "outs/sco-{sim},{sco}.rds"
    log:    "logs/sco-{sim},{sco}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# calculate feature selection
rule calc_sel:
    priority: 96
    input:  "code/03-sel.R",
            "code/03-sel-{sel}.R",
            x = res_sco_by_sim
    params: lambda wc, input: ";".join(input.x)
    output: "outs/sel-{sim},{sel}.rds"
    log:    "logs/sel-{sim},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {params} {output[0]}" {input[0]} {log}'''

# reprocessing using selected features
rule rep_data:
    priority: 95
    input:  "code/04-rep.R",
            "data/01-fil/{sim}.rds",
            rules.calc_sel.output
    output: "data/02-rep/{sim},{sel}.rds"
    log:    "logs/rep_data-{sim},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
            {input[2]} {output}" {input[0]} {log}'''

# calculate clustering evaluation statistics
rule calc_sta:
    priority: 94
    input:  "code/05-sta.R",
            "code/05-sta-{sta}.R",
            rules.rep_data.output
    output: "outs/sta-{sim},{sel},{sta}.rds"
    log:    "logs/sta-{sim},{sel},{sta}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# differential abundance/analysis 
rule run_das:
    priority: 94
    input:  "code/06-das.R",
            "code/06-das-{das}.R",
            rules.rep_data.output
    output: "outs/das-{sim},{sel},{das}.rds"
    log:    "logs/das-{sim},{sel},{das}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# calculate performance of detect true markers

rule calc_eva:
    priority: 95
    input:  "code/07-eva.R",
            "code/07-eva-roc.R",
            rules.calc_sco.output
    output: "outs/eva-{sim},{sco},{eva}.rds"
    log:    "logs/eva-{sim},{sco},{eva}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''





# VISUALIZATION ========================================================

for val in VAL:
    rule:
        priority: 49
        input:  expand("code/08-plot_{val}-{{plt}}.R", val = val), x = res[val]
        params: lambda wc, input: ";".join(input.x)
        output: expand("plts/{val}-{{plt}}.pdf", val = val)
        log:    expand("logs/plot_{val}-{{plt}}.Rout", val = val)
        shell:  '''
            {R} CMD BATCH --no-restore --no-save "--args\
            {params} {output[0]}" {input[0]} {log}'''

#rule plot_roc:
#    priority: 49
#    input:  "code/08-plot-roc_curve.R", x = res_roc
#    params: lambda wc, input: ";".join(input.x)
#    output: "plts/roc.pdf"
#    log:    "logs/plot-roc.Rout"
#    shell:  '''
#        {R} CMD BATCH --no-restore --no-save "--args\
#        {params} {output[0]}" {input[0]} {log}'''
