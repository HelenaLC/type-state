import json
import itertools

configfile: "config.yaml"
R = config["R"]

SCO = glob_wildcards("code/02-sco-type_{x}.R").x
SEL = glob_wildcards("code/03-sel-{x}.R").x
STA = glob_wildcards("code/05-sta-{x}.R").x

# magnitude of type, state & batch effect
T = list(range(0, 120, 20))
S = list(range(0, 120, 20))
B = [0]

SIM = ["t{},s{},b{}".format(t,s,b) for t in T for s in S for b in B]

res_sim = expand(
    "data/00-sim/t{t},s{s},b{b}.rds", 
    t = T, s = S, b = B)
res_fil = expand(
    "data/01-fil/{sim}.rds", 
    sim = SIM)
res_sco = expand(
    "outs/sco-{sim},{sco}.rds",
    sim = SIM, sco = SCO)
res_sel = expand(
    "outs/sel-{sim},{sco},{sel}.rds",
    sim = SIM, sco = SCO, sel = SEL)
res_rep = expand(
    "data/02-rep/{sim},{sco},{sel}.rds", 
    sim = SIM, sco = SCO, sel = SEL)
res_sta = expand(
    "outs/sta-{sim},{sco},{sel},{sta}.rds",
    sim = SIM, sco = SCO, sel = SEL, sta = STA)
res_roc = expand(
    "outs/roc-{sim},{sco},{sel}.rds",
    sim = SIM, sco = SCO, sel = SEL)

# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
        "session_info.txt",
        # simulation, pre- & re-processing
        res_sim, res_fil, res_rep,
        # scoring & evaluation
        res_sco, res_sel, res_sta, res_roc,
        # visualization
        expand(
            "plts/{plt}.pdf", 
            plt = ["pca", "sco", "sta", "roc_curve"])

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
            "code/02-sco-type_{sco}.R",
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
           rules.calc_sco.output
   output: "outs/sel-{sim},{sco},{sel}.rds"
   log:    "logs/sel-{sim},{sco},{sel}.Rout"
   shell: '''
       {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
       {input[1]} {input[2]} {output[0]}" {input[0]} {log}'''

# reprocessing using selected features
rule rep_data:
    priority: 94
    input:  "code/04-rep.R",
            rules.fil_data.output,
            rules.calc_sel.output
    output: "data/02-rep/{sim},{sco},{sel}.rds"
    log:    "logs/rep_data-{sim},{sco},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# calculate clustering evaluation statistics
rule calc_sta:
    priority: 95
    input:  "code/05-sta.R",
            "code/05-sta-{sta}.R",
            rules.rep_data.output
    output: "outs/sta-{sim},{sco},{sel},{sta}.rds"
    log:    "logs/sta-{sim},{sco},{sel},{sta}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

# calculate performance of detect true markers

rule det_tm:
    priority: 95
    input:  "code/06-roc.R",
            "code/06-roc-type.R",
            rules.calc_sel.output
    output: "outs/roc-{sim},{sco},{sel}.rds"
    log:    "logs/roc-{sim},{sco},{sel}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''

# COLLECTION ===========================================================

def res_sta_by_sco(wildcards):
    return expand("outs/sta-{sim},{sco},{sel},{sta}.rds",
        sim = SIM, sco = wildcards.sco, sel = SEL, sta = STA)

# VISUALIZATION ========================================================

rule plot_pca:
    priority: 49
    input:  "code/08-plot_pca.R", 
            x = expand("data/01-fil/{sim}-cd.rds", sim = SIM)
    params: lambda wc, input: ";".join(input.x)
    output: "plts/pca.pdf"
    log:    "logs/plot-pca.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''

rule plot_sco:
    priority: 49
    input:  "code/08-plot-sco.R", x = res_sco
    params: lambda wc, input: ";".join(input.x)
    output: "plts/sco.pdf"
    log:    "logs/plot-sco.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''

rule plot_sta:
    priority: 49
    input:  "code/08-plot-sta.R", x = res_sta
    params: lambda wc, input: ";".join(input.x)
    output: "plts/sta.pdf"
    log:    "logs/plot-sta.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''

rule plot_roccurve:
    priority: 49
    input:  "code/08-plot-roccurve.R", x = res_roc
    params: lambda wc, input: ";".join(input.x)
    output: "plts/roc_curve.pdf"
    log:    "logs/plot-roccurve.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''

