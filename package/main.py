"""
Author: Marios Nikolaidis
Git: marios-nikolaidis
email: marionik23@gmail.com
"""

import argparse
import config
from multiprocessing import cpu_count

parser = argparse.ArgumentParser(
        prog="package",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="This is a description of the program",
        epilog="This is a description of the program",
        allow_abbrev=False
)
subparsers = parser.add_subparsers(help="The various modules of the program")

#install = parser.add_argument_group("install")
#install.add_argument(
#        "--install",
#        action="store_true",
#        help="Install the package",
#)

species_demarcation = subparsers.add_parser("species_demarcation", help="species_demarcation module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
species_demarcation.add_argument("-in", metavar="in" , help="Genome fasta directory", required=True)
species_demarcation.add_argument("-o", metavar="out", help="Output directory", required=True)
species_demarcation_fastani = species_demarcation.add_argument_group("FastANI options")
species_demarcation_fastani.add_argument("--fastani_cores", help="Number of cores for FastANI", default=config.species_demarcation_cores)
species_demarcation_fastani.add_argument("--kmer", help="kmer size", default=config.species_demarcation_kmer)
species_demarcation_fastani.add_argument("--fragLen", help="Fragment length", default=config.species_demarcation_fraglen)
species_demarcation_fastani.add_argument("--minFraction", help="Minimum fraction", default=config.species_demarcation_minfraction)
species_demarcation_mcl = species_demarcation.add_argument_group("FastANI options")
species_demarcation_mcl.add_argument("--inflation", help="Inflation parameter for MCL", default=config.species_demarcation_mcl_inflation)
species_demarcation_mcl.add_argument("--mcl_cores", help="Number of cores for MCL", default=config.species_demarcation_cores)
#
orthologues = subparsers.add_parser("orthologues", help="orthologues module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
eggnog.add_argument("-in", help="Input core proteins file", required=True)
eggnog.add_argument("-fasta_dir", help="Input fasta directory", required=True)
eggnog.add_argument("-o", help="Output directory")
eggnog.add_argument("--eggnog_results", help="Pre-computed eggnog results")
eggnog_mapper = eggnog.add_argument_group("eggNOG mapper options")
eggnog_mapper.add_argument("--cores", help="Number of cores", default=config.eggnog_cores, type=int)
eggnog_mapper.add_argument("--pident", help="Percent identity", default=config.eggnog_pident)
eggnog_mapper.add_argument("--qcov", help="Query coverage", default=config.eggnog_qcov)
eggnog_mapper.add_argument("--scov", help="Suject coverage", default=config.eggnog_scov)
eggnog_installer = eggnog.add_argument_group("Installer options")
eggnog.add_argument("--install", help="Install eggnog-mapper and the databases", action="store_true")

# smbgc module
smbgc = subparsers.add_parser("smbgc", help="smbgc module", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
smbgc.add_argument("-fasta_dir", metavar="fasta_dir", help="Input fasta directory")
smbgc.add_argument("-o", metavar="output_dir", help="Output directory")
smbgc.add_argument("--cores", metavar="N", help="Number of cores", default=config.smbgc_cores, type=int)
smbgc.add_argument("--strictness", metavar="strictness", help="antiSMASH strictness settings", default=config.smbgc_strictness, type=str)
smbgc.add_argument("--genefinding_tool", metavar="tool", help="Gene finding tool used by antiSMASH", default=config.smbgc_genefinding_tool, type=str)
#smbgc.add_argument("--remote", help="Submit queries to antiSMASH web service", action="store_true")

args = parser.parse_args()
print(args)
