import datetime
import csv
import multiprocessing
import os
import logging
import argparse
from Pipeline import Instance
from Phylogenomics import Phylogenomics_instance
import sys
from Bio import SeqIO
import re
import pathlib

parser = argparse.ArgumentParser(usage='python main.py PARAM_FILE', epilog=' Wrapper for the application')
parser.add_argument("parameters_file", metavar='f', help=" Parameters file location")
args = parser.parse_args()

param_file = args.parameters_file

def _defaultParams():
	system = sys.platform
	input_dir = os.getcwd()
	output_dir = ""
	threads_to_use =  multiprocessing.cpu_count() - 2

	params = {
		"sys": system,
		"in": pathlib.Path(input_dir),
		"out": pathlib.Path(output_dir),
		"fasta_files" : [],
		"evalue": 1e-5, 
		"max_target_seqs": 1,
		"num_threads": threads_to_use,
		"SD": 2,
		"bootstrap": 500,
		"distance": "Kimura"
	}

	return params

def _updateParams(param_file, params):

	param_csv = csv.reader(open(param_file, "r"), delimiter="\t")
	for lines in param_csv:
		key, value = lines
		params[key] = value
	
	if params["out"] == "": params["out"] = params["in"]
	params["in"] = pathlib.Path(params["in"])
	params["out"] = pathlib.Path(params["out"])
	# Gather the fasta files for the analysis
	fasta_files = list((params["in"] / "Fasta_files").glob("*"))
	params["fasta_files"] = fasta_files
	
	return params



def _readEvolfile(*args, **kwargs):
	"""
	Read the tab delimited volutionary groups file
	Org \t Evol group
	"""
	evol_groups = {}
	evol_group_stats = {}
	evol_groups_f = kwargs['in'] / "evol_groups.txt"
	if os.path.exists(evol_groups_f) == True:
		fin = open(evol_groups_f)
		for lines in csv.reader(fin, delimiter = "\t" ):
			org, evogroup = lines
			evol_groups[org] = evogroup
			if evogroup not in evol_group_stats:
				evol_group_stats[evogroup] = 1
			else:
				evol_group_stats[evogroup] += 1

	return evol_groups, evol_group_stats

def _readInstancefile(*args, **kwargs):
	"""
	Instance name\tReference
	"""
	instances = {}
	instances_f = kwargs['in'] / "instances_file.txt"
	with open(instances_f, "r") as fin:
		for lines in csv.reader(fin, delimiter = "\t"):
			inst_name, ref = lines
			instances[inst_name] = ref
	return instances


if __name__ == "__main__":
	now = datetime.datetime.now()
	start_dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
	print("Script started at ", start_dt_string)

	default_params = _defaultParams()
	execution_params = _updateParams(param_file, default_params)
	execution_params['out'].mkdir(exist_ok=True)
	logfile = execution_params['out'] / "LOGFILE.log"
	logging.basicConfig(filename=logfile, format='%(asctime)s %(message)s', level=logging.DEBUG)
	evol_groups, evol_group_stats = _readEvolfile(**execution_params)
	instances = _readInstancefile(**execution_params)
	# print(execution_params['in'])
	full_instance_counter = 1
	for instance_name in instances:
		ref = instances[instance_name]

		now = datetime.datetime.now()
		dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
		print( "####Starting Pipeline for Instance %s at %s #####" %(instance_name, dt_string))
		logging.info(" Initilizing instance " + str(instance_name))
		instance_n = Instance(execution_params, evol_groups, evol_group_stats, instance_name, ref)
		#instance_n.create_dirs()
		#instance_n.reciprocal_blast()
		#instance_n.parse_blast_results()
		instance_n.create_cog_matrix()
		
		instance_n_phylo = Phylogenomics_instance(execution_params, ref, instance_name)
		#instance_n_phylo.createCOGFasta()
		#instance_n_phylo.alignCOGs()
		instance_n_phylo.createSupersequence()
		instance_n_phylo.filterAln() 
		instance_n_phylo.computeTree(iqtree=True) # Doesn't work properly in WSL

end_dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print("The scripts started at: %s \n Finished at %s" %(start_dt_string, end_dt_str))
