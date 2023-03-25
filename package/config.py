# Configuration file for the pipeline
# By M. Nikolaidis

from sys import platform as platform
from multiprocessing import cpu_count

system = platform 

## 1. Species demarcation configs
species_demarcation_cores = cpu_count() - 2
species_demarcation_kmer = 16
species_demarcation_fraglen = 3000
species_demarcation_minfraction = 0.2
species_demarcation_mcl_inflation = 2

## 2. Orthologues configs
orthologues_cores = 6
orthologues_evalue = 1e-5
orthologues_SD_filter_value = 2

#def _updateParams(param_file, params):
#
#	for lines in param_csv:
#        param_csv = csv.reader(open(param_file, "r"), delimiter="\t")
#		key, value = lines
#		params[key] = value
#	
#	if params["out"] == "": params["out"] = params["in"]
#	params["in"] = pathlib.Path(params["in"])
#	params["out"] = pathlib.Path(params["out"])
#	# Gather the fasta files for the analysis
#	fasta_files = list((params["in"] / "Fasta_files").glob("*"))
#	params["fasta_files"] = fasta_files
#	
#	return params
#
