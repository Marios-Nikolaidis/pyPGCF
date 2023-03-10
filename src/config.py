# Configuration file for the pipeline
# By M. Nikolaidis

from sys import platform as platform
from pathlib import Path
from multiprocessing import cpu_count

system = platform 
BLAST_EVALUE = 1e-5
BLAST_THREADS = 6
PHYLOGENOMIC_THREADS = cpu_count() - 2
SD = 2

def _updateParams(param_file, params):

	for lines in param_csv:
        param_csv = csv.reader(open(param_file, "r"), delimiter="\t")
		key, value = lines
		params[key] = value
	
	if params["out"] == "": params["out"] = params["in"]
	params["in"] = pathlib.Path(params["in"])
	params["out"] = pathlib.Path(params["out"])
	# Gather the fasta files for the analysis
	fasta_files = list((params["in"] / "Fasta_files").glob("*"))
	params["fasta_files"] = fasta_files
	
	return params

