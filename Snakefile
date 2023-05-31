import json
import itertools

configfile: "config.yaml"
R = config["R"]

# get datasets
DATSETS = glob_wildcards("code/scripts/00-get_data-{x}.R").x
SCORES = glob_wildcards("code/scripts/02-score-{x}.R").x
N_FEATURES = glob_wildcards("code/scripts/03-sel_{x}.R").x

#TYPE_SCORES = glob_wildcards("output/score/").x
#DATATYPES = glob_wildcards("code/scripts/01-fil_data-{x}.R").x

# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule all: 
    input:
        expand([
        # Load data, preprocessing, initial clustering
            "data/00-raw/{datset}.rds",
            "data/01-fil/{datset}.rds",
        # stateness and typeness score
            "output/score/{datset},{score}.rds",
        
        # Feature selection with new scores
            "data/02-new/{datset},{score},{n_features}_new.rds",
            "output/sel/{datset},{score},{n_features}_sel.rds",

        # UMAP with newly selected features
            "output/UMAP/{datset},{score},{n_features}_UMAP.png",

        # evaluation scores
            "output/eval/{datset},{score},{n_features}_eval.rds"],
            datset = DATSETS,
            score = SCORES,
            n_features = N_FEATURES),
        

rule get_data:
    priority: 98
	input: 	"code/scripts/00-get_data.R",
			"code/scripts/00-get_data-{datset}.R"
	output:	"data/00-raw/{datset}.rds"
	log:	"logs/get_data-{datset}.Rout"
	shell:	'''
        {R} CMD BATCH --no-restore --no-save\
        "--args {input[1]} {output[0]}" {input[0]} {log}'''

rule fil_data:
    priority: 97
    input:  "code/scripts/01-fil_data_splatter.R",
                    rules.get_data.output
    output: "data/01-fil/{datset}.rds"
    log:    "logs/fil_data-{datset}.Rout"
    shell:  '''
        {R} CMD BATCH --no-restore --no-save\
            "--args {input[1]} {output[0]}" {input[0]} {log}'''


rule cal_score:
    priority: 96
    input:  "code/scripts/02-score.R",
            "code/scripts/02-score-{score}.R",
            rules.fil_data.output
    output: "output/score/{datset},{score}.rds"
    log:    "logs/score-{datset},{score}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {output}" {input[0]} {log}'''

rule new_data:
    priority: 95
    input:  "code/scripts/03-sel.R",
            "code/scripts/03-sel_{n_features}.R",
            rules.fil_data.output,
            rules.cal_score.output
    output: "data/02-new/{datset},{score},{n_features}_new.rds",
            "output/sel/{datset},{score},{n_features}_sel.rds"
    log:    "logs/new-{datset},{score},{n_features}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {input[2]} {input[3]} {output[0]} {output[1]}" {input[0]} {log}'''

rule eval_scores:
    priority: 94
    input: "code/scripts/04-eval.R",
            rules.new_data.output
    output: "output/eval/{datset},{score},{n_features}_eval.rds"
    log:    "logs/eval-{datset},{score},{n_features}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {output}" {input[0]} {log}'''

rule plot_UMAPs:
    priority: 93
    input: "code/scripts/05-plot-UMAPs.R",
            rules.new_data.output
    output: "output/UMAP/{datset},{score},{n_features}_UMAP.png"
    log:    "logs/UMAP-{datset},{score},{n_features}.Rout"
    shell: '''
        {R} CMD BATCH --no-restore --no-save "--args {input[1]}\
        {output}" {input[0]} {log}'''

