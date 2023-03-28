import pandas as pd
from numpy import nan as np_nan
from pathlib import Path
from typing import List
import os
import logging
#import re
#import multiprocessing
#from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def _execute_cmd(cmd: str):
    os.system(cmd)

class PhylogenomicInstance():
    def __init__(self, orthology_matrix_f: Path, cores: int, system: str, fasta_dir: Path, out_dir: Path, resources_dir: Path):
        orthology_matrix = pd.read_csv(orthology_matrix_f, sep="\t", index_col=0)
        self.system = system
        self.fasta_dir = fasta_dir
        self.out_dir = out_dir / "phylogenomic_tree"
        self.og_fasta_dir = self.out_dir / "OGs_fasta"
        self.og_fasta_dir_aln = self.out_dir / "OGs_fasta_aln"
        self.iqtree_results_dir = self.out_dir / "iqtree"
        self.orthology_matrix = orthology_matrix
        self.ref = orthology_matrix.index.name
        self.genomes = orthology_matrix.columns.tolist()
        self.cores = cores
        self.muscle_bin_path = resources_dir / "muscle3.8.31"
        self.gblocks_bin_path = resources_dir / "Gblocks"
        self.iqtree_bin_path = resources_dir / "iqtree2"

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)
        self.og_fasta_dir.mkdir(exist_ok=True, parents=True)
        self.og_fasta_dir_aln.mkdir(exist_ok=True, parents=True)
        self.iqtree_results_dir.mkdir(exist_ok=True, parents=True)

    def create_og_fasta(self):
        """
        Create a fasta file for each cluster of orthologous genes
        :return: None
        """
        self.orthology_matrix.applymap(lambda x: np_nan if x == "X" else x)
        cog_matrix = self.orthology_matrix.dropna()
        organisms: List = [self.ref]
        organisms.extend(cog_matrix.columns.tolist())
        fasta_files = self.fasta_dir.glob("*")
        fasta_files_dict = {f.stem:f for f in fasta_files}
        for org_count, organism in enumerate(organisms):
            fasta_file = fasta_files_dict[organism]
            protein_records = SeqIO.to_dict(SeqIO.parse(open(fasta_file,"r"), format="fasta"))
            if org_count == 0:
                print(fasta_file)
                genes = list(cog_matrix.index)
            else:
                genes = list(cog_matrix[organism].values)
            for x, gene in enumerate(genes):
                info = protein_records[gene]
                info.description = ""
                info.name = ""
                info.id = info.id + "-" + organism
                fout = self.og_fasta_dir / ("OG" + str(x) + ".fa")
                fout_handle = open(fout, "a")
                SeqIO.write(info, fout_handle, "fasta")
                fout_handle.close()
    # Write presequity for dash in gene name. Need to find a way to fix this

    def align_og_fasta(self):
        """
        Use muscle to align each file of orthologous group
        """
        orthologous_groups_files = list(self.og_fasta_dir.glob("*"))
        commands = []
        for file in orthologous_groups_files:
            cmd = " ".join([str(self.muscle_bin_path),
                            "-in",
                            str(file),
                            "-out",
                            str(self.og_fasta_dir_aln / file.name),
                            "-quiet"
                        ])
            commands.append(cmd)
        with ProcessPoolExecutor(self.cores) as executor:
            _ = executor.map(_execute_cmd,commands)
    
    def create_supersequence_file(self):
        """
        Join the aligned COGs into a supersequence
        :return: None
        """
        # def _sortAlphanum(iteratable):
        #     # Helper func Sort the given motif list alphanumerically :return: sorted list 
        #     int_convert = lambda text: int(text) if text.isdigit() else text 
        #     sorting_key = lambda key: [ int_convert(c) for c in re.split('([0-9]+)', key) ] 
        #     return sorted(iteratable, key = sorting_key)
        
        def init_supersequence_file(ref, genomes, superseq_file):
            """
            Initialise the supersequence file with empty sequences for each organism
            """
            seqrecord_list = []
            organisms = [ref]
            organisms.extend(genomes[:])
            for org in organisms:
                record = SeqRecord(Seq(""), id = org, name = org, description="")
                seqrecord_list.append(record)
            superseq_file_handle_out = open(superseq_file, "w")
            SeqIO.write(seqrecord_list,superseq_file_handle_out, "fasta")
    
        superseq_file = self.out_dir / "supersequence.fa"
        init_supersequence_file(self.ref, self.genomes, superseq_file)
        superseq_records = SeqIO.to_dict(AlignIO.read(str(superseq_file), "fasta"))

        aln_files = self.fasta_dir.glob("*")
        # aln_file = [str(f) for f in aln_files]
        # aln_files = _sortAlphanum(aln_files) 
        # To know that this is the correct order of genes
        for aln_file in aln_files:
            aln_records = SeqIO.to_dict(AlignIO.read(str(aln_file), "fasta"))
            aln_records = {records.split("-")[1] : aln_records[records] for records in aln_records} # Rename the keys, need to use the org name
            aln_organisms = list(aln_records.keys())
            for organism in superseq_records:
                if organism in aln_organisms:
                    superseq_records[organism].seq = superseq_records[organism].seq + aln_records[organism].seq
        superseq_seqrecords = list(superseq_records.values())
        superseq_file_handle_out = open(superseq_file, "w")
        SeqIO.write(superseq_seqrecords, superseq_file_handle_out, "fasta")
    
    
    def filter_aln(self):
        """ Filter the supersequence alignment using Gblocks with default parameters """
        cmd = " ".join([str(self.gblocks_bin_path), 
                        str(self.out_dir / "supersequence.fa"),
                        "-s=y -e=-gb -p=y" 
                        ])
        _execute_cmd(cmd)
    
    def compute_tree(self):
        """
        Compute the phylogenomic tree using IQtree2
        """
        supersequence_file = self.out_dir / "supersequence.fa-gb"
        cmd = "".join([str(self.iqtree_bin_path),
                        " -m TEST -merit AIC -alrt 1000 -T ",
                        str(self.cores), " -s ", str(supersequence_file)])
        _execute_cmd(cmd)

    def run_phylogenomic(self) -> int:
        self.setup_directories()
        self.create_og_fasta()
        logging.debug("Created a file for each orthologous group. Dir %s"  %(self.out_dir))
        logging.debug(" Initiating alignment of individual orthologous groups ")
        self.align_og_fasta()
        logging.debug(" Finished COG alignments in dir %s" %(self.out_dir))
        self.create_supersequence_file()
        logging.debug("Created the supersequence file")
        self.filter_aln()
        logging.debug(" Filtered the alignment with GBlocks")
        logging.debug(f" Initiating computation of phylogenomic tree")
        self.compute_tree()
        logging.debug("Computed Phylogenomic Tree ")
        return 0
