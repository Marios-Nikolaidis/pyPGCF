#!/usr/bin/env python3
"""
Author: Marios Nikolaidis
Git: marios-nikolaidis
email: marionik23@gmail.com
"""

import argparse
import config
from pathlib import Path
from species_demarcation.species_demarcation import SpeciesDemarcator
from orthologues.orthologues import Orthologues_identifier
from core.core import Core_identifier
from phylogenomic.phylogenomic import Phylogenomic
from eggNOG.eggNOG import eggNOGRunner, eggNOGParser
from smbgc.smbgc import smBGCLocalRunner, smBGCParser, smBGCInstaller

parser = argparse.ArgumentParser(
        prog="package",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This is a description of the program",
        epilog="This is a description of the program",
        allow_abbrev=False
)
subparsers = parser.add_subparsers(help="The various modules of the program", dest="module", required=True)

species_demarcation = subparsers.add_parser("species_demarcation", help="species_demarcation module", formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Assign the input genomes to species clusters using FastANI and MCL")
species_demarcation.add_argument("-in", metavar="input" , help="Genome fasta directory", required=True)
species_demarcation.add_argument("-o", metavar="out", help="Output directory", required=True)
species_demarcation_fastani = species_demarcation.add_argument_group("FastANI options")
species_demarcation_fastani.add_argument("--fastani_cores", help="Number of cores for FastANI", default=config.species_demarcation_cores)
species_demarcation_fastani.add_argument("--kmer", help="kmer size", default=config.species_demarcation_kmer)
species_demarcation_fastani.add_argument("--fraglen", help="Fragment length", default=config.species_demarcation_fraglen)
species_demarcation_fastani.add_argument("--minfraction", help="Minimum fraction", default=config.species_demarcation_minfraction)
species_demarcation_mcl = species_demarcation.add_argument_group("MCL options")
species_demarcation_mcl.add_argument("--inflation", help="Inflation parameter for MCL", default=config.species_demarcation_mcl_inflation)
species_demarcation_mcl.add_argument("--mcl_cores", help="Number of cores for MCL", default=config.species_demarcation_cores)
#
orthologues = subparsers.add_parser("orthologues", help="orthologues module", formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Calculate the orthologues of the dataset using one reference strain")
orthologues_basic = orthologues.add_argument_group("Basic options")
orthologues_basic.add_argument("-in", metavar="in", help="Fasta directory", required=True)
orthologues_basic.add_argument("-o", metavar="out", help="Output directory", required=True)
ref_exclusive = orthologues_basic.add_mutually_exclusive_group(required=True)
ref_exclusive.add_argument("-ref", help="Reference strain")
ref_exclusive.add_argument("-ref_list", help="List of reference strains to use")
orthologues_basic.add_argument("--type", help="Type of input [prot (DIAMOND) or nucl(BLASTN)]", default="prot")
orthologues_blast = orthologues.add_argument_group("DIAMOND/BLASTN options")
orthologues_blast.add_argument("--cores", help="Number of cores", default=config.orthologues_cores)
orthologues_blast.add_argument("--evalue", help="E-value cut-off", default=config.orthologues_evalue)
orthologues_blast.add_argument("--dmnd_sensitivity", help="Sensitivity settings for DIAMOND", default=config.orthologues_dmnd_sensitivity)
orthologues_blast.add_argument("--no_filter_orthologues", help="Do not filter orthologues", action="store_true")

core = subparsers.add_parser("core", help="core module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
core.add_argument("-in", help="Input orthology matrix", required=True)
core.add_argument("-o", help="Output directory", required=True)
core.add_argument("--species", help="Input species assignment matrix (excel format)")
core.add_argument("--core_perc", help="Percent presence of a gene in a cluster to be considered core", default=config.core_core_perc)

phylogenomic = subparsers.add_parser("phylogenomic", help="phylogenomic module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
phylogenomic.add_argument("-fasta_dir", help="Input fasta directory", required=True)
phylogenomic.add_argument("-in", help="Input Orthology matrix", required=True)
phylogenomic.add_argument("-o", help="Output directory", required=True)
phylogenomic.add_argument("--cores", help="Number of cores", default=config.phylogenomic_cores)
phylogenomic_alns = phylogenomic.add_argument_group("Intermediate fasta files options")
phylogenomic_alns.add_argument("--no_keep_fasta", help="Remove the orthologous groups fasta files after completion", action="store_true")
phylogenomic_tree = phylogenomic.add_argument_group("IQTree2 options")
phylogenomic_tree.add_argument("--tree_model", help="Tree model")

# eggnog module
eggnog = subparsers.add_parser("eggnog", help="eggnog module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
eggnog.add_argument("-in", metavar="in", help="Input core proteins file", required=True)
eggnog.add_argument("-fasta_dir", metavar="fasta_dir", help="Input fasta directory", required=True)
eggnog.add_argument("-o", metavar="o", help="Output directory", required=True)
eggnog.add_argument("--eggnog_results", metavar="file", help="Pre-computed eggnog results")
eggnog_mapper = eggnog.add_argument_group("eggNOG mapper options")
eggnog_mapper.add_argument("--cores", metavar="N", help="Number of cores", default=config.eggnog_cores, type=int)
eggnog_mapper.add_argument("--pident", metavar="N", help="Percent identity", default=config.eggnog_pident)
eggnog_mapper.add_argument("--qcov", metavar="N", help="Query coverage", default=config.eggnog_qcov)
eggnog_mapper.add_argument("--scov", metavar="N", help="Suject coverage", default=config.eggnog_scov)
eggnog_installer = eggnog.add_argument_group("Installer options")
eggnog.add_argument("--install", help="Install eggnog-mapper and the databases", action="store_true")

# smbgc module
smbgc = subparsers.add_parser("smbgc", help="smbgc module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
smbgc.add_argument("-fasta_dir", metavar="fasta_dir", help="Input fasta directory")
smbgc.add_argument("-o", metavar="output_dir", help="Output directory")
smbgc.add_argument("--cores", metavar="N", help="Number of cores", default=config.smbgc_cores, type=int)
smbgc.add_argument("--strictness", metavar="strictness", help="antiSMASH strictness settings", default=config.smbgc_strictness, type=str)
smbgc.add_argument("--genefinding_tool", metavar="tool", help="Gene finding tool used by antiSMASH", default=config.smbgc_genefinding_tool, type=str)
smbgc.add_argument("--install", help="Install antiSMASH and the required database", action="store_true")
#smbgc.add_argument("--remote", help="Submit queries to antiSMASH web service", action="store_true")

args = vars(parser.parse_args())
if args["module"] == "species_demarcation":
    in_dir = Path(args["in"])
    out_dir = Path(args["o"])
    fastani_cores = args["fastani_cores"]
    kmer = args["kmer"]
    fraglen = args["fraglen"]
    minfraction = args["minfraction"]
    inflation = args["inflation"]
    mcl_cores = args["mcl_cores"]
    demarcator = SpeciesDemarcator(in_dir, out_dir, fastani_cores, kmer, fraglen, minfraction, inflation, mcl_cores)
    demarcator.assign_species() # Need to fix this

if args["module"] == "orthologues":
    fasta_in_dir = Path(args["in"])
    out_dir = Path(args["o"])
    ref = args["ref"]
    ref_list = args["ref_list"]
    input_type = args["type"]
    cores = args["cores"]
    evalue = args["evalue"]
    dmdnd_sensitivity = args["dmnd_sensitivity"]
    no_filter_orthologues = args["no_filter_orthologues"]
    orthologues_identifier = Orthologues_identifier(fasta_in_dir, out_dir, ref, ref_list, input_type, cores, evalue, dmdnd_sensitivity, no_filter_orthologues)
    orthologues_identifier.calculate_orthologues()

if args["module"] == "core":
    og_matrix_in = Path(args["in"])
    out_dir = Path(args["o"])
    species_file = args["species"]
    if species_file != None:
        species_file = Path(species_file)
    core_perc = args["core_perc"]
    core_identifier = Core_identifier(og_matrix_in, out_dir, species_file, core_perc)
    core_identifier.calculate_core()

if args["module"] == "phylogenomic":
    fasta_dir = Path(args["fasta_dir"])
    og_matrix_in = Path(args["in"])
    out_dir = Path(args["o"])
    cores = args["cores"]
    no_keep_fasta = args["no_keep_fasta"]
    tree_model = args["tree_model"]
    phylogenomic = Phylogenomic(og_matrix_in, cores, fasta_dir, out_dir, no_keep_fasta, tree_model)
    phylogenomic.run_phylogenomic()

if args["module"] == "eggnog":
    fasta_dir = Path(args["fasta_dir"])
    core_proteins_file = Path(args["in"])
    out_dir = Path(args["o"])
    precomputed_results = args["eggnog_results"]
    if precomputed_results != None:
        precomputed_results = Path(precomputed_results)
    cores = args["cores"]
    pident = args["pident"]
    qcov = args["qcov"] # Query coverage
    scov = args["scov"] # Subject coverage
    install = args["install"]
    if install:
        pass
    else:
        runner = eggNOGRunner(fasta_dir, core_proteins_file, out_dir, cores, pident, qcov, scov)
        runner.execute_eggnog_mapper()
        parser = eggNOGParser(fasta_dir, core_proteins_file, out_dir)
        parser.gather_eggnog_results()

if args["module"] == "smbgc":
    install = args["install"]
    if install:
        installer = smBGCInstaller()
        installer.install_antismash()
        installer.install_databases()
    else:
        fasta_dir = Path(args["fasta_dir"])
        out_dir = Path(args["o"])
        cores = args["cores"]
        strictness = args["strictness"]
        genefinding_tool = args["genefinding_tool"]
        local_runner = smBGCLocalRunner(fasta_dir, out_dir, cores, strictness, genefinding_tool)
        local_runner.analyze_genomes()
        parser = smBGCParser(out_dir, cores)
        parser.gather_results()
