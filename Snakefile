import json
import itertools

configfile: "config.yaml"
R = config["R"]

<<<<<<< HEAD
# get datasets
TYPE_SCORES = glob_wildcards("code/scripts/02-score-type{x}.R").x
STATE_SCORES = glob_wildcards("code/scripts/02-score-state{x}.R").x
SEL_TYPE = glob_wildcards("code/scripts/03-selType_{x}.R").x
TYPE = ["70", "75", "80", "85", "90", "95", "100"]
STATE = ["0", "5", "10", "15", "20", "25", "30"]
B = ["0"]

=======
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
>>>>>>> origin/main

# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
<<<<<<< HEAD
        expand([
        # Load data, preprocessing, initial clustering
            "data/00-sim/sim-t{t}-s{s}-b{b}.rds",
            "data/01-fil/fil-t{t}-s{s}-b{b}.rds",
            "outs/sim_pca-t{t}-s{s}-b{b}.png",

        # stateness and typeness score
            "outs/type{t_score}-t{t}-s{s}-b{b}.rds",
            "outs/state{s_score}-t{t}-s{s}-b{b}.rds"],
        
        # Feature selection with new scores
        #    "data/02-new/new-t{t}-s{s}-b{b},{t_score},type_{sel_type}.rds",
        #    "outs/sel-t{t}-s{s}-b{b},{t_score},type_{sel_type}.rds"],

        # UMAP with newly selected features
        #    "output/UMAP/{datset},{score},{n_features}_UMAP.png",

        # evaluation scores
        #    "output/eval/{datset},{score},{n_features}_eval.rds"],
        #    datset = DATSETS,
        #    score = SCORES,
        #    n_features = N_FEATURES,
            t = TYPE,
            s = STATE,
            b = B,
            t_score = TYPE_SCORES,
            s_score = STATE_SCORES
            #sel_type = SEL_TYPE
            )
        
=======
        "session_info.txt",
        # simulation, pre- & re-processing
        res_sim, res_fil, res_rep,
        # scoring & evaluation
        res_sco, res_sel, res_sta,
        # visualization
        expand(
            "plts/{plt}.pdf", 
            plt = ["pca", "sco", "sta"])

rule session_info:
    priority: 100
    input:  "code/10-session_info.R"
    output: "session_info.txt"
    log:    "logs/session_info.Rout"
    shell:  '''
    {R} CMD BATCH --no-restore --no-save\
    "--args {output}" {input} {log}'''
>>>>>>> origin/main

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

# calculate evaluation statistics
rule calc_sta:
    priority: 95
    input:  "code/05-sta.R",
            "code/05-sta-{sta}.R",
            rules.rep_data.output
    output: "outs/sta-{sim},{sco},{sel},{sta}.rds"
    log:    "logs/sta-{sim},{sco},{sel},{sta}.Rout"
    shell: '''
<<<<<<< HEAD
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''

rule state_score:
    priority: 95
    input:  "code/scripts/02-score.R",
            "code/scripts/02-score-state{s_score}.R",
            rules.fil_data.output
    output: "outs/state{s_score}-t{t}-s{s}-b{b}.rds"
    log:    "logs/state{s_score}-t{t}-s{s}-b{b}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''   

rule sole_sel_type:
    priority: 95
    input:  "code/scripts/03-selType.R",
            "code/scripts/03-selType_{sel_type}.R",
            rules.fil_data.output,
            rules.type_score.output
    output: "data/02-new/new-t{t}-s{s}-b{b},{t_score},type_{sel_type}.rds",
            "outs/sel-t{t}-s{s}-b{b},{t_score},type_{sel_type}.rds"
    log:    "logs/new-t{t}-s{s}-b{b},{t_score},type{sel_type}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args fun={input[1]}\
        sce={input[2]} score={input[3]} new={output[0]} sel={output[1]}" {input[0]} {log}'''

#rule eval_scores:
#    priority: 94
#    input: "code/scripts/04-eval.R",
#            rules.new_data.output
#    output: "output/eval/{datset},{score},{n_features}_eval.rds"
#    log:    "logs/eval-{datset},{score},{n_features}.Rout"
#    shell: '''
#        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
#        {output}" {input[0]} {log}'''

#rule plot_UMAPs:
#    priority: 93
#    input: "code/scripts/05-plot-UMAPs.R",
#            rules.new_data.output
#    output: "output/UMAP/{datset},{score},{n_features}_UMAP.png"
#    log:    "logs/UMAP-{datset},{score},{n_features}.Rout"
#    shell: '''
#        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
#       {output}" {input[0]} {log}'''

# PLOTS ================================================================
rule plot_sim_pca:
    priority: 96
    input: "code/scripts/01-sim-pca.R",
            rules.fil_data.output
    output: "outs/sim_pca-t{t}-s{s}-b{b}.png"
    log:   "logs/sim_pca-t{t}-s{s}-b{b}.Rout"
=======
        {R} CMD BATCH --no-restore --no-save "--args wcs={wildcards}\
        {input[1]} {input[2]} {output}" {input[0]} {log}'''

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
>>>>>>> origin/main
    shell:  '''
        {R} CMD BATCH --no-restore --no-save "--args\
        {params} {output[0]}" {input[0]} {log}'''