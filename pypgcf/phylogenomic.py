from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from os import system
from pathlib import Path
from typing import Generator, List, Union

from numpy import nan as np_nan
from tqdm import tqdm

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd


def _execute_cmd(cmd: str):
    system(cmd)


class Phylogenomic:
    def __init__(
        self,
        orthology_matrix_f: Path,
        cores: int,
        fasta_dir: Path,
        out_dir: Path,
        no_keep_fasta: bool,
        tree_model: Union[str, None],
        debug: bool = False,
    ):
        orthology_matrix = pd.read_csv(orthology_matrix_f, sep="\t", index_col=0)
        self.fasta_dir = fasta_dir
        self.out_dir = out_dir / "Phylogenomic_tree"
        self.og_fasta_dir = self.out_dir / "OGs_fasta"
        self.og_fasta_dir_aln = self.out_dir / "OGs_fasta_aln"
        self.iqtree_results_dir = self.out_dir / "iqtree"
        self.orthology_matrix = orthology_matrix
        self.ref = orthology_matrix.index.name
        self.genomes = orthology_matrix.columns.tolist()
        self.cores = cores
        self.no_keep_fasta = no_keep_fasta
        self.debug = debug
        self.tree_model = tree_model
        if tree_model == None:
            self.tree_model = "TEST"

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)
        self.og_fasta_dir.mkdir(exist_ok=True, parents=True)
        self.og_fasta_dir_aln.mkdir(exist_ok=True, parents=True)
        self.iqtree_results_dir.mkdir(exist_ok=True, parents=True)

    def _replace_empty_with_na_in_orthology_matrix(self):
        self.orthology_matrix = self.orthology_matrix.applymap(
            lambda x: np_nan if x == "X" else x
        )

    def create_og_fasta(self):
        """
        Create a fasta file for each cluster of orthologous genes
        :return: None
        """
        self._replace_empty_with_na_in_orthology_matrix()
        self.orthology_matrix = self.orthology_matrix.dropna()
        if self.orthology_matrix.shape[0] == 0:
            print(
                f"Could not identify orthologous groups with orthologues from all genomes"
            )
            import sys

            sys.exit()
        organisms: List = [self.orthology_matrix.index.name]
        organisms.extend(self.orthology_matrix.columns.tolist())
        fasta_files: Generator = self.fasta_dir.glob("*")
        fasta_files_dict = {f.stem: f for f in fasta_files}
        for org_count, organism in enumerate(organisms):
            fasta_file = fasta_files_dict[organism]
            protein_records = SeqIO.to_dict(
                SeqIO.parse(open(fasta_file, "r"), format="fasta")
            )
            if org_count == 0:
                genes = list(self.orthology_matrix.index)
            else:
                genes = list(self.orthology_matrix[organism].values)
            for x, gene in enumerate(genes):
                info = protein_records[gene]
                info.description = organism
                info.name = ""
                info.id = info.id
                fout = self.og_fasta_dir / ("OG" + str(x) + ".fa")
                fout_handle = open(fout, "a")
                SeqIO.write(info, fout_handle, "fasta")
                fout_handle.close()

    def align_og_fasta(self):
        """
        Use muscle to align each file of orthologous group
        """
        orthologous_groups_files = list(self.og_fasta_dir.glob("*"))
        commands = []
        for file in orthologous_groups_files:
            cmd = " ".join(
                [
                    "muscle",
                    "-in",
                    str(file),
                    "-out",
                    str(self.og_fasta_dir_aln / file.name),
                    "-quiet",
                ]
            )
            commands.append(cmd)
        with ProcessPoolExecutor(self.cores) as executor:
            list(
                tqdm(
                    executor.map(_execute_cmd, commands),
                    total=len(commands),
                    desc="Aligning OG fasta files",
                    ascii=True,
                    leave=True,
                )
            )

    def create_superalignment_file(self):
        """
        Join the aligned orthologous groups into a superalignment
        :return: None
        """

        def init_superalignment_file(ref, genomes, superseq_file):
            """
            Initialise the superalignment file with empty sequences for each organism
            """
            seqrecord_list = []
            organisms = [ref]
            organisms.extend(genomes[:])
            for org in organisms:
                record = SeqRecord(Seq(""), id=org, name=org, description="")
                seqrecord_list.append(record)
            superseq_file_handle_out = open(superseq_file, "w")
            SeqIO.write(seqrecord_list, superseq_file_handle_out, "fasta")

        superseq_file = self.out_dir / "superalignment.fa"
        init_superalignment_file(self.ref, self.genomes, superseq_file)
        superseq_records = SeqIO.to_dict(AlignIO.read(str(superseq_file), "fasta"))
        aln_files = self.og_fasta_dir_aln.glob("*")

        for aln_file in tqdm(
            aln_files, desc="Joining orthologous groups", ascii=True, leave=True
        ):
            # aln_records = SeqIO.to_dict(AlignIO.read(str(aln_file), "fasta"))
            # aln_records = {records.description : aln_records[records] for records in aln_records} # Rename the keys, need to use the org name
            parser = SeqIO.parse(str(aln_file), "fasta")
            aln_records = {
                record.description.replace(record.name + " ", ""): record
                for record in parser
            }  # Rename the keys, need to use the org name
            aln_organisms = list(aln_records.keys())
            for organism in superseq_records:
                if organism in aln_organisms:
                    superseq_records[organism].seq += aln_records[organism].seq
        superseq_seqrecords = list(superseq_records.values())
        superseq_file_handle_out = open(superseq_file, "w")
        SeqIO.write(superseq_seqrecords, superseq_file_handle_out, "fasta")

    def filter_superalignment(self):
        """Filter the superalignment using Gblocks with default parameters"""
        cmd = " ".join(
            ["Gblocks", str(self.out_dir / "superalignment.fa"), "-s=y -e=-gb -p=y"]
        )
        _execute_cmd(cmd)
        htm_file = self.out_dir / "superalignment.fa-gb.htm"
        htm_file.unlink()

    def compute_tree(self):
        """
        Compute the phylogenomic tree using IQtree2
        """
        superalignment_file = self.out_dir / "superalignment.fa-gb"
        cmd = " ".join(
            [
                "iqtree2",
                f"-m {self.tree_model}",
                "--quiet",
                "-merit AIC",
                "-alrt 1000",
                f"-T {str(self.cores)}",
                f"-s {superalignment_file}",
            ]
        )
        _execute_cmd(cmd)

    def move_iqtree_files(self):
        files = list(self.out_dir.glob("superalignment.fa-gb.*"))
        for f in files:
            renamed = self.iqtree_results_dir / f.name
            if f.suffix == ".treefile":
                renamed = self.iqtree_results_dir / "superalignment_IQTree2.nwk"
            f.rename(renamed)

    def clean_fasta_files(self):
        directories_to_clean = [self.og_fasta_dir, self.og_fasta_dir_aln]
        for directory in directories_to_clean:
            files = directory.glob("*")
            for f in files:
                f.unlink()
            directory.unlink()

    def run_phylogenomic(self) -> Union[None, int]:
        self.setup_directories()
        print(
            f"Creating fasta files of each orthologous group: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.create_og_fasta()
        print(
            f"Aligning fasta files of each orthologous group: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.align_og_fasta()
        print(
            f"Creating super-alignment file: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.create_superalignment_file()
        print(
            f"Filtering super-alignment file: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.filter_superalignment()
        print(
            f"Computing phylogenomic tree: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}"
        )
        self.compute_tree()
        self.move_iqtree_files()
        print(f"Done: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        if self.no_keep_fasta:
            self.clean_fasta_files()
        if self.debug:
            return 0
