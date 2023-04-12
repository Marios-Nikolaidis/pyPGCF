"""
Author: Marios Nikolaidis
Git: https://github.com/Marios-Nikolaidis
email: marionik23@gmail.com
"""
import argparse
from . import config
from pathlib import Path
from .species_demarcation import SpeciesDemarcator
from .orthologues import Orthologues_identifier
from .core import Core_identifier
from .phylogenomic import Phylogenomic
from .eggnog import eggNOGRunner, eggNOGParser, eggNOGInstaller
from .smbgc import smBGCLocalRunner, smBGCParser, smBGCInstaller


def main():
    parser = argparse.ArgumentParser(
        prog="package",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="pyPGCF: PhyloGenomic, Core and Fingerprint analysis software",
        allow_abbrev=False,
    )
    subparsers = parser.add_subparsers(
        help="The various modules of the program", dest="module", required=True
    )

    species_demarcation = subparsers.add_parser(
        "species_demarcation",
        help="species_demarcation module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Assign the input genomes to species clusters using FastANI and MCL",
    )
    species_demarcation.add_argument(
        "-in", metavar="input", help="Genome fasta directory", required=True
    )
    species_demarcation.add_argument(
        "-o", metavar="out", help="Output directory", required=True
    )
    species_demarcation_fastani = species_demarcation.add_argument_group(
        "FastANI options"
    )
    species_demarcation_fastani.add_argument(
        "--fastani_cores",
        help="Number of cores for FastANI",
        default=config.species_demarcation_cores,
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--kmer",
        metavar="N",
        help="kmer size (<= 16)",
        default=config.species_demarcation_kmer,
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--fraglen",
        metavar="N",
        help="Fragment length",
        default=config.species_demarcation_fraglen,
        type=int,
    )
    species_demarcation_fastani.add_argument(
        "--minfraction",
        metavar="N",
        help="Minimum fraction",
        default=config.species_demarcation_minfraction,
        type=float,
    )
    species_demarcation_mcl = species_demarcation.add_argument_group("MCL options")
    species_demarcation_mcl.add_argument(
        "--inflation",
        help="Inflation parameter for MCL",
        default=config.species_demarcation_mcl_inflation,
    )
    species_demarcation_mcl.add_argument(
        "--mcl_cores",
        help="Number of cores for MCL",
        default=config.species_demarcation_cores,
        type=int,
    )

    # orthologues module
    orthologues = subparsers.add_parser(
        "orthologues",
        help="orthologues module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Calculate the orthologues of the dataset using one reference strain",
    )
    orthologues_basic = orthologues.add_argument_group("Basic options")
    orthologues_basic.add_argument(
        "-in", metavar="in", help="Fasta directory", required=True
    )
    orthologues_basic.add_argument(
        "-o", metavar="out", help="Output directory", required=True
    )
    ref_exclusive = orthologues_basic.add_mutually_exclusive_group(required=True)
    ref_exclusive.add_argument("-ref", help="Reference strain")
    ref_exclusive.add_argument("-ref_list", help="List of reference strains to use")
    orthologues_basic.add_argument(
        "--type", help="Type of input [prot (DIAMOND) or nucl(BLASTN)]", default="prot"
    )
    orthologues_blast = orthologues.add_argument_group("DIAMOND/BLASTN options")
    orthologues_blast.add_argument(
        "--cores", help="Number of cores", default=config.orthologues_cores
    )
    orthologues_blast.add_argument(
        "--evalue", help="E-value cut-off", default=config.orthologues_evalue
    )
    orthologues_blast.add_argument(
        "--dmnd_sensitivity",
        help="Sensitivity settings for DIAMOND",
        default=config.orthologues_dmnd_sensitivity,
    )
    orthologues_blast.add_argument(
        "--no_filter_orthologues", help="Do not filter orthologues", action="store_true"
    )

    # core module
    core = subparsers.add_parser(
        "core",
        help="core module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    core_input_group = core.add_mutually_exclusive_group(required=True)
    core_input_group.add_argument("-in", help="Input orthology matrix")
    core_input_group.add_argument(
        "-ref_list",
        help="List of orthology matrices and species annotation files to use for multiple analyses",
    )
    core.add_argument("-o", help="Output directory", required=True)
    core.add_argument(
        "--species", help="Input species assignment matrix (excel format)"
    )
    core.add_argument(
        "--core_perc",
        help="Percent presence of a protein/gene in a cluster to be considered core",
        default=config.core_core_perc,
    )

    # phylogenomic module
    phylogenomic = subparsers.add_parser(
        "phylogenomic",
        help="phylogenomic module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    phylogenomic.add_argument("-fasta_dir", help="Input fasta directory", required=True)
    phylogenomic.add_argument("-in", help="Input Orthology matrix", required=True)
    phylogenomic.add_argument("-o", help="Output directory", required=True)
    phylogenomic.add_argument(
        "--cores", help="Number of cores", default=config.phylogenomic_cores
    )
    phylogenomic_alns = phylogenomic.add_argument_group(
        "Intermediate fasta files options"
    )
    phylogenomic_alns.add_argument(
        "--no_keep_fasta",
        help="Remove the orthologous groups fasta files after completion",
        action="store_true",
    )
    phylogenomic_tree = phylogenomic.add_argument_group("IQTree2 options")
    phylogenomic_tree.add_argument(
        "--tree_model", help="Specific evolutionary model for tree calculation"
    )

    # eggnog module
    eggnog = subparsers.add_parser(
        "eggnog",
        help="eggnog module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    eggnog_input = eggnog.add_mutually_exclusive_group()
    eggnog_input.add_argument(
        "-in", metavar="in", help="Excel file with core proteins (from 'core' module)"
    )
    eggnog_input.add_argument(
        "-in_list",
        metavar="in",
        help="File with a list of excel files (for multiple reference strains)",
    )
    eggnog.add_argument("-fasta_dir", metavar="fasta_dir", help="Input fasta directory")
    eggnog.add_argument("-o", metavar="o", help="Output directory")
    # eggnog.add_argument("--eggnog_results", metavar="file", help="Pre-computed eggnog results")
    eggnog_mapper = eggnog.add_argument_group("eggNOG mapper options")
    eggnog_mapper.add_argument(
        "--cores",
        metavar="N",
        help="Number of cores",
        default=config.eggnog_cores,
        type=int,
    )
    eggnog_mapper.add_argument(
        "--pident", metavar="N", help="Percent identity", default=config.eggnog_pident
    )
    eggnog_mapper.add_argument(
        "--qcov", metavar="N", help="Query coverage", default=config.eggnog_qcov
    )
    eggnog_mapper.add_argument(
        "--scov", metavar="N", help="Suject coverage", default=config.eggnog_scov
    )
    eggnog_installer = eggnog.add_argument_group("Installer options")
    eggnog_installer.add_argument(
        "--install", help="Install the eggNOG database", action="store_true"
    )

    # smbgc module
    smbgc = subparsers.add_parser(
        "smbgc",
        help="smbgc module",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    smbgc.add_argument("-fasta_dir", metavar="fasta_dir", help="Input fasta directory")
    smbgc.add_argument("-o", metavar="output_dir", help="Output directory")
    smbgc.add_argument(
        "--cores",
        metavar="N",
        help="Number of cores",
        default=config.smbgc_cores,
        type=int,
    )
    smbgc.add_argument(
        "--strictness",
        metavar="strictness",
        help="antiSMASH strictness settings",
        default=config.smbgc_strictness,
        type=str,
    )
    smbgc.add_argument(
        "--genefinding_tool",
        metavar="tool",
        help="Gene finding tool used by antiSMASH",
        default=config.smbgc_genefinding_tool,
        type=str,
    )
    smbgc.add_argument(
        "--install",
        help="Install antiSMASH and the required database",
        action="store_true",
    )
    # smbgc.add_argument("--remote", help="Submit queries to antiSMASH web service", action="store_true")

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
        demarcator = SpeciesDemarcator(
            in_dir,
            out_dir,
            fastani_cores,
            kmer,
            fraglen,
            minfraction,
            inflation,
            mcl_cores,
        )
        demarcator.assign_species()  # Need to fix this

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
        orthologues_identifier = Orthologues_identifier(
            fasta_in_dir,
            out_dir,
            ref,
            ref_list,
            input_type,
            cores,
            evalue,
            dmdnd_sensitivity,
            no_filter_orthologues,
        )
        orthologues_identifier.calculate_orthologues()

    if args["module"] == "core":
        out_dir = Path(args["o"])
        species_file = args["species"]
        if species_file != None:
            species_file = Path(species_file)
        og_matrix_in = args["in"]
        og_matrix_list_f = args["ref_list"]
        if og_matrix_list_f == None:
            og_matrix_list = [og_matrix_in]
        else:
            og_matrix_list = []
            with open(og_matrix_list_f, "r") as rf:
                for line in rf:
                    line = line.rstrip()
                    og_matrix_list.append(line)
        og_matrix_list = [Path(og_matrix_in) for og_matrix_in in og_matrix_list]
        for og_matrix_in in og_matrix_list:
            core_perc = args["core_perc"]
            core_identifier = Core_identifier(
                og_matrix_in, out_dir, species_file, core_perc
            )
            core_identifier.calculate_core()

    if args["module"] == "phylogenomic":
        fasta_dir = Path(args["fasta_dir"])
        og_matrix_in = Path(args["in"])
        out_dir = Path(args["o"])
        cores = args["cores"]
        no_keep_fasta = args["no_keep_fasta"]
        tree_model = args["tree_model"]
        phylogenomic = Phylogenomic(
            og_matrix_in, cores, fasta_dir, out_dir, no_keep_fasta, tree_model
        )
        phylogenomic.run_phylogenomic()

    if args["module"] == "eggnog":
        install = args["install"]
        if install:
            installer = eggNOGInstaller()
            installer.download_databases()
            return
        fasta_dir = args["fasta_dir"]
        out_dir = args["o"]
        if fasta_dir == None:
            print(f"-fasta_dir is needed")
            return
        if out_dir == None:
            print(f"-o is needed")
            return
        fasta_dir = Path(fasta_dir)
        out_dir = Path(out_dir)
        cores = args["cores"]
        pident = args["pident"]
        qcov = args["qcov"]  # Query coverage
        scov = args["scov"]  # Subject coverage
        core_proteins_file = args["in"]
        core_protein_files_reflist = args["in_list"]
        if core_proteins_file == None and core_protein_files_reflist == None:
            print("-in or -in_list are needed")
            return
        if core_proteins_file != None:
            core_proteins_file_list = [core_proteins_file]
        else:
            core_proteins_file_list = []
            with open(core_protein_files_reflist, "r") as rf:
                for line in rf:
                    core_proteins_file_list.append(line.rstrip())
        core_proteins_file_list = [
            Path(core_proteins_file) for core_proteins_file in core_proteins_file_list
        ]
        for core_proteins_file in core_proteins_file_list:
            core_proteins_file = Path(core_proteins_file)
            runner = eggNOGRunner(
                fasta_dir, core_proteins_file, out_dir, cores, pident, qcov, scov
            )
            runner.execute_eggnog_mapper()
        parser = eggNOGParser(fasta_dir, core_proteins_file_list, out_dir)
        parser.gather_eggnog_results()

    if args["module"] == "smbgc":
        install = args["install"]
        if install:
            installer = smBGCInstaller()
            installer.install_databases()
        else:
            fasta_dir = Path(args["fasta_dir"])
            out_dir = Path(args["o"])
            cores = args["cores"]
            strictness = args["strictness"]
            genefinding_tool = args["genefinding_tool"]
            local_runner = smBGCLocalRunner(
                fasta_dir, out_dir, cores, strictness, genefinding_tool
            )
            local_runner.analyze_genomes()
            parser = smBGCParser(out_dir, cores)
            parser.gather_results()
