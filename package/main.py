"""
Author: Marios Nikolaidis
Git: marios-nikolaidis
email: marionik23@gmail.com
"""

import argparse
import config

parser = argparse.ArgumentParser(
        prog="package",
        description="This is a description of the program",
        epilog="This is a description of the program",
        allow_abbrev=False
)
subparsers = parser.add_subparsers(dest="K", help="The various modules of the program")

#install = parser.add_argument_group("install")
#install.add_argument(
#        "--install",
#        action="store_true",
#        help="Install the package",
#)

species_demarcation = subparsers.add_parser("species_demarcation", help="species demarcation")
species_demarcation.add_argument("--in", help="Genome fasta directory")
species_demarcation.add_argument("--o", help="Output directory")
species_demarcation.add_argument("--cores", help="Number of cores", default=config.species_demarcation_cores)
species_demarcation.add_argument("--kmer", help="kmer size for fastANI", default=config.species_demarcation_kmer)
species_demarcation.add_argument("--fragLen", help="Fragment length for fastANI", default=config.species_demarcation_fraglen)
species_demarcation.add_argument("--minFraction", help="Minimum fraction for fastANI", default=config.species_demarcation_minfraction)
species_demarcation.add_argument("--inflation", help="Inflation parameter for mcl", default=config.species_demarcation_mcl_inflation)
#
orthologues = subparsers.add_parser("orthologues", help="orthologues")
orthologues_basic = orthologues.add_argument_group("Basic options")
orthologues_basic.add_argument("--i", help="Fasta directory")
orthologues_basic.add_argument("--o", help="Output directory")
orthologues_basic.add_argument("--type", help="Type of input", default="protein")
orthologues_basic.add_argument("--ref", help="Reference strain")
orthologues_blast = orthologues.add_argument_group("BLAST/DIAMOND options")
orthologues_blast.add_argument("--cores", help="Number of cores", default=config.orthologues_cores)
orthologues_blast.add_argument("--evalue", help="E-value cut-off", default=config.orthologues_evalue)
# TODO: Complete

core = subparsers.add_parser("core", help="core")
core.add_argument("--cog_in", help="Input COG matrix")
core.add_argument("--o", help="Output directory")
core.add_argument("--clust_in", help="Input species-cluster matrix", keep=True)
core.add_argument("--core_perc", help="Percent presence of a gene in a cluster to be considered core", default=100)

phylogenomic = subparsers.add_parser("phylogenomic", help="phylogenomic")
phylogenomic.add_argument("--fasta_dir", help="Input fasta directory")
phylogenomic.add_argument("--cog_in", help="Input COG matrix")
phylogenomic.add_argument("--o", help="Output directory")
phylogenomic.add_argument("--cores", help="Number of cores", default=cpu_count()-2)
phylogenomic_tree = phylogenomic.add_argument_group("IQTree2 options")
phylogenomic_tree.add_argument("--tree_model", help="Tree model")

# eggnog module
eggnog = subparsers.add_parser("eggnog", help="eggnog")
eggnog.add_argument("--install", help="Install eggnog-mapper", action="store_true") # Not yet implemented
eggnog.add_argument("--ref_fasta", help="Input reference fasta")
eggnog.add_argument("--o", help="Output directory")
eggnog.add_argument("--core", help="Core proteins table")
eggnog.add_argument("--eggnog_results", help="Pre-computed eggnog results", action="store_true")
eggnog.add_argument("--cores", help="Number of cores", default=cpu_count()-2)

# smbgc module
smbgc = subparsers.add_parser("smbgc", help="smbgc")
smbgc.add_argument("--i", help="Input fasta directory")
smbgc.add_argument("--o", help="Output directory")
smbgc.add_argument("--cores", help="Number of cores", default=cpu_count()-2)
smbgc.add_argument("--remote", help="Submit queries to antiSMASH web service", action="store_true")


args = parser.parse_args()

