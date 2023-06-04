import json
import itertools

configfile: "config.yaml"
R = config["R"]

# get datasets
#DATSETS = glob_wildcards("code/scripts/00-get_data-{x}.R").x
TYPE_SCORES = glob_wildcards("code/scripts/02-score-type{x}.R").x
STATE_SCORES = glob_wildcards("code/scripts/02-score-state{x}.R").x
N_FEATURES = glob_wildcards("code/scripts/03-sel_{x}.R").x
TYPE = ["0", "25", "50", "75", "100"]
STATE = ["0", "25", "50", "75", "100"]
BCV = ["1"]

#TYPE_SCORES = glob_wildcards("output/score/").x
#DATATYPES = glob_wildcards("code/scripts/01-fil_data-{x}.R").x

# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
        expand([
        # Load data, preprocessing, initial clustering
            "data/00-sim/sim-t{t}-s{s}-b{b}.rds",
            "data/01-fil/fil-t{t}-s{s}-b{b}.rds",
            "output/plots/sim_pca/pca-t{t}-s{s}-b{b}.png",

        # stateness and typeness score
            "output/score/data-t{t}-s{s}-b{b},type{t_score}.rds",
            "output/score/data-t{t}-s{s}-b{b},state{s_score}.rds"],
        
        # Feature selection with new scores
        #    "data/02-new/{datset},{score},{n_features}_new.rds",
        #    "output/sel/{datset},{score},{n_features}_sel.rds",

        # UMAP with newly selected features
        #    "output/UMAP/{datset},{score},{n_features}_UMAP.png",

        # evaluation scores
        #    "output/eval/{datset},{score},{n_features}_eval.rds"],
        #    datset = DATSETS,
        #    score = SCORES,
        #    n_features = N_FEATURES,
            t = TYPE,
            s = STATE,
            b = BCV,
            t_score = TYPE_SCORES,
            s_score = STATE_SCORES
            )
        

rule sim_data:
    priority: 98
    input: "code/scripts/00-sim_data.R",
    output: "data/00-sim/sim-t{t}-s{s}-b{b}.rds"
    log:	"logs/sim-t{t}-s{s}-b{b}.Rout"
    params: t = lambda wildcards: wildcards.t,
            s = lambda wildcards: wildcards.s,
            b = lambda wildcards: wildcards.b
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args t={params.t}\
        s={params.s} b={params.b} res={output}" {input[0]} {log}'''


rule fil_data:
    priority: 97
    input:  "code/scripts/01-fil_data_splatter.R",
            rules.sim_data.output
    output: "data/01-fil/fil-t{t}-s{s}-b{b}.rds"
    log:    "logs/fil_data-t{t}-s{s}-b{b}.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save\
            "--args {input[1]} {output[0]}" {input[0]} {log}'''


rule type_score:
    priority: 95
    input:  "code/scripts/02-score.R",
            "code/scripts/02-score-type{t_score}.R",
            rules.fil_data.output
    output: "output/score/data-t{t}-s{s}-b{b},type{t_score}.rds"
    log:    "logs/score-data-t{t}-s{s}-b{b},type{t_score}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''

rule state_score:
    priority: 95
    input:  "code/scripts/02-score.R",
            "code/scripts/02-score-state{s_score}.R",
            rules.fil_data.output
    output: "output/score/data-t{t}-s{s}-b{b},state{s_score}.rds"
    log:    "logs/score-data-t{t}-s{s}-b{b},state{s_score}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''        

#rule new_data:
#    priority: 95
#    input:  "code/scripts/03-sel.R",
#            "code/scripts/03-sel_{n_features}.R",
#            rules.fil_data.output,
#            rules.type_score.output
#    output: "data/02-new/data-t{t}-s{s}-b{b},{score},{n_features}_new.rds",
#            "output/sel/data-t{t}-s{s}-b{b},{score},{n_features}_sel.rds"
#    log:    "logs/new-data-t{t}-s{s}-b{b},{score},{n_features}.Rout"
#    shell: '''
#        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
#        {input[2]} {input[3]} {output[0]} {output[1]}" {input[0]} {log}'''

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
    output: "output/plots/sim_pca/pca-t{t}-s{s}-b{b}.png"
    log:   "logs/sim_pca-t{t}-s{s}-b{b}.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save\
            "--args {input[1]} {output[0]}" {input[0]} {log}'''