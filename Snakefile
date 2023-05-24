import json
import itertools


configfile: "config.yaml"
R = config["R"]

# get datasets, subsets & methods
DATSETS = glob_wildcards("code/scripts/00-get_data-{x}.R").x

# PREPROCESSING ================================================================

# reproducibly retrieve dataset from public source
rule get_data:
    priority: 98
	input: 	"code/scripts/00-get_data.R",
			"code/scripts/00-get_data-{datset}.R"
	output:	"data/00-raw/{datset}.rds"
	log:	"logs/get_data-{datset}.Rout"
	shell:	'''
	{R} CMD BATCH --no-restore --no-save\
