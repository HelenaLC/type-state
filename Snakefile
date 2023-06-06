import json
import itertools

configfile: "config.yaml"
R = config["R"]

# get datasets
TYPE_SCORES = glob_wildcards("code/scripts/02-score-type{x}.R").x
STATE_SCORES = glob_wildcards("code/scripts/02-score-state{x}.R").x
SEL_TYPE = glob_wildcards("code/scripts/03-selType_{x}.R").x
TYPE = ["70", "75", "80", "85", "90", "95", "100"]
STATE = ["0", "5", "10", "15", "20", "25", "30"]
B = ["0"]


# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
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
        

rule sim_data:
    priority: 98
    input: "code/scripts/00-sim_data.R",
    output: "data/00-sim/sim-t{t}-s{s}-b{b}.rds"
    log:	"logs/sim-t{t}-s{s}-b{b}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save\
        "--args wcs={wildcards} res={output}" {input[0]} {log}'''


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
    output: "outs/type{t_score}-t{t}-s{s}-b{b}.rds"
    log:    "logs/type{t_score}-t{t}-s{s}-b{b}.Rout"
    shell: '''
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
    shell:  '''
        {R} CMD BATCH --no-restore --no-save\
            "--args {input[1]} {output[0]}" {input[0]} {log}'''