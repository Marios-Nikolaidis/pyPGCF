"""Add functional annotation to the proteins in the core genes excels"""
import pandas as pd
from collections import Counter
from Bio import SeqIO
from pathlib import Path
import os

class eggNOGInstaller():
    def __init__(self, resources_dir: Path, debug: bool):
        self.resources_dir = resources_dir
        self.debug = debug

    def install_emapper(self):
        print("Installing eggnog mapper")
        cmd = "pip install eggnog-mapper"
        if self.debug:
            cmd = "pip install -v eggnog-mapper"
        os.system(cmd)
    
    def download_databases(self):
        cmd = f"download_eggnog_data.py -y --data_dir {self.resources_dir}"
        if self.debug:
            cmd = f"download_eggnog_data.py -s -y --data_dir {self.resources_dir}"
        os.system(cmd)

    def install(self):
        self.install_emapper()
        self.download_databases()

class eggNOGRunner():
    def __init__(self, fasta_dir: Path, core_protein_table_f: Path, out_dir: Path, cores: int, pident: float, qcov: float, scov: float):
        self.fasta_dir = fasta_dir
        self.core_protein_table = pd.read_excel(core_protein_table_f, index_col=0)
        self.out_dir = out_dir / "eggNOG"
        self.cores = cores
        self.pident = pident
        self.qcov =  qcov
        self.scov = scov

    def get_reference_fasta_file(self):
        ref = self.core_protein_table.index.name
        found = False
        for fasta_file in self.fasta_dir.glob("*"):
            if ref not in fasta_file.name:
                continue
            found = True
            self.fasta_file = fasta_file
        if not found:
            raise Exception(f"{ref} fasta file not found in {self.fasta_dir.name}")

    def create_eggnog_cmd(self):
        eggnog_cmd = " ".join([
            "emapper.py",
            f"--cpu {self.cores}",
            f"--pident {self.pident}",
            f"--query_cov {self.qcov}",
            f"--subject_cov {self.scov}",
            f"-o {self.out_dir/'eggNOG_results.csv'}",
            f"-i {self.fasta_file}"
            ])
        self.eggnog_cmd = eggnog_cmd

    def execute_eggnog_mapper(self):
        self.get_reference_fasta_file()
        self.create_eggnog_cmd()
        os.system(self.eggnog_cmd)

class eggNOGParser():
    def __init__(self, fasta_dir: Path, core_protein_table_f: Path, out_dir: Path):
        self.fasta_dir = fasta_dir
        self.core_protein_table = pd.read_excel(core_protein_table_f, index_col=0)
        self.out_dir = out_dir / "eggNOG"
        self.eggnog_headers = ["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", 
            "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
            "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
            "CAZy", "BiGG_Reaction", "PFAMs"
        ]
        self.category_annotation = {
            "A":"RNA processing and modification",
            "B":"Chromatin structure and dynamics",
            "C":"Energy production and conversion",
            "D":"Cell cycle control, cell division, chromosome partitioning",
            "E":"Amino acid transport and metabolism",
            "F":"Nucleotide transport and metabolism",
            "G":"Carbohydrate transport and metabolism",
            "H":"Coenzyme transport and metabolism",
            "I":"Lipid transport and metabolism",
            "J":"Translation, ribosomal structure and biogenesis",
            "K":"Transcription",
            "L":"Replication, recombination and repair",
            "M":"Cell wall/membrane/envelope biogenesis",
            "N":"Cell motility",
            "O":"Post-translational modification, protein turnover, and chaperones",
            "P":"Inorganic ion transport and metabolism",
            "Q":"Secondary metabolites biosynthesis, transport, and catabolism",
            "S":"Function unknown",
            "T":"Signal transduction mechanisms",
            "U":"Intracellular trafficking, secretion, and vesicular transport",
            "V":"Defense mechanisms",
            "W":"Extracellular structures",
            "Z":"Cytoskeleton",
        }

    def get_reference_fasta_file(self):
        ref = self.core_protein_table.index.name
        found = False
        for fasta_file in self.fasta_dir.glob("*"):
            if ref not in fasta_file.name:
                continue
            found = True
            self.fasta_file = fasta_file
        if not found:
            raise Exception(f"{ref} fasta file not found in {self.fasta_dir.name}")

    def load_all_reference_proteins(self) -> dict:
        record = SeqIO.parse(str(self.fasta_file), "fasta")
        return {protein.id: ["S"] for protein in record}
    
    def get_category(self, x: str) -> list:
        if x == "-":
            return ["S"]
        return list(x)
    
    def load_eggnog_results(self) -> dict:
        eggnog_results_f = self.out_dir / "eggNOG_results.csv"
        annotation_df = pd.read_csv(eggnog_results_f, sep="\t", index_col=0, comment="#", names=self.eggnog_headers)
        annotation_df["COG_category"] = annotation_df["COG_category"].apply(
            lambda x: self.get_category(x)
        )
        return annotation_df["COG_category"].to_dict()

    def parse_eggnog_results(self):
        # Load reference proteins
        ref_data = self.load_all_reference_proteins()
        # Load eggnog annotation
        eggnog_annotation = self.load_eggnog_results()
        # Update reference data 
        ref_data.update(eggnog_annotation)
        self.proteome_cog = ref_data

    def get_cog_annotation_of_protein_subset(self):
        pass
         
    def hypergeometric_test(self):
        pass
        
        #fout = "Genome_functional_categories_20221117.xlsx"
        #data_df = pd.DataFrame.from_dict(data_dictionary, orient="index")
        #data_df.to_excel(fout)
