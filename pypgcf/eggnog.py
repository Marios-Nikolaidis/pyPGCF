"""Add functional annotation to the proteins in the core genes excels"""
import pandas as pd
from scipy.stats import hypergeom
from datetime import datetime
from Bio import SeqIO
from pathlib import Path
from typing import Tuple, List
import os

class eggNOGInstaller():
    def __init__(self, debug: bool=False):
        self.debug = debug

    # def install_emapper(self):
    #     print("Installing eggnog mapper")
    #     cmd = "pip install eggnog-mapper"
    #     if self.debug:
    #         cmd = "pip install --dry-run -v eggnog-mapper"
    #     os.system(cmd)
    
    def download_databases(self):
        cmd = f"download_eggnog_data.py -y"
        #cmd = f"download_eggnog_data.py -y --data_dir {self.resources_dir}"
        if self.debug:
            #cmd = f"download_eggnog_data.py -s -y --data_dir {self.resources_dir}"
            cmd = f"download_eggnog_data.py -s -y"
        os.system(cmd)

class eggNOGRunner():
    def __init__(self, fasta_dir: Path, core_protein_table_f: Path, out_dir: Path, cores: int, pident: float, qcov: float, scov: float, debug: bool=False):
        #self.resources_dir = resources_dir
        self.fasta_dir = fasta_dir
        self.core_protein_table = pd.read_excel(core_protein_table_f, index_col=0)
        self.ref = str(self.core_protein_table.index.name)
        self.out_dir = out_dir / "eggNOG"
        self.eggnog_raw_results_file = self.out_dir / (self.ref + "_eggNOG_results.csv")
        self.cores = cores
        self.pident = pident
        self.qcov =  qcov
        self.scov = scov
        self.debug = debug

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)

    def get_reference_fasta_file(self):
        found = False
        for fasta_file in self.fasta_dir.glob("*"):
            if self.ref != fasta_file.stem:
                continue
            found = True
            self.fasta_file = fasta_file
        if not found:
            raise FileNotFoundError(f"{self.ref} fasta file not found in {self.fasta_dir.name}")

    def create_eggnog_cmd(self):
        eggnog_cmd = " ".join([
            "emapper.py",
            #f"--data_dir {self.resources_dir}",
            f"--cpu {self.cores}",
            f"--pident {self.pident}",
            f"--query_cov {self.qcov}",
            f"--subject_cov {self.scov}",
            f"-o {self.eggnog_raw_results_file}",
            f"-i {self.fasta_file}",
            " > /dev/null "
            ])
        self.eggnog_cmd = eggnog_cmd

    def clean_unessecary_output(self):
        files_to_remove = [
            self.eggnog_raw_results_file.with_suffix(".csv.emapper.hits"),
            self.eggnog_raw_results_file.with_suffix(".csv.emapper.seed_orthologs"),
        ]
        file_to_rename = self.eggnog_raw_results_file.with_suffix(".csv.emapper.annotations")
        for f in files_to_remove:
            f.unlink()
        file_to_rename.rename(self.eggnog_raw_results_file)

    def execute_eggnog_mapper(self):
        print(f"Executing eggNOG-mapper: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}, reference strain: {self.ref}")
        self.setup_directories()
        self.get_reference_fasta_file()
        self.create_eggnog_cmd()
        self.setup_directories()
        if self.debug:
            status = os.system(self.eggnog_cmd)
            self.execute_status =status
        else:
            os.system(self.eggnog_cmd)
        self.clean_unessecary_output()

class eggNOGParser():
    def __init__(self, fasta_dir: Path, core_protein_table_file_list: List[Path], out_dir: Path):
        self.fasta_dir = fasta_dir
        self.core_protein_table_file_list = core_protein_table_file_list
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
        self.ref = str(self.core_protein_table.index.name)
        found = False
        for fasta_file in self.fasta_dir.glob("*"):
            if self.ref not in fasta_file.name:
                continue
            found = True
            self.fasta_file = fasta_file
        if not found:
            raise Exception(f"{self.ref} fasta file not found in {self.fasta_dir.name}")

    def load_all_reference_proteins(self) -> dict:
        record = SeqIO.parse(str(self.fasta_file), "fasta")
        return {protein.id: ["S"] for protein in record}
    
    def turn_category_to_list(self, x: str) -> list:
        if x == "-":
            return ["S"]
        return list(x)
    
    def get_assigned_cog_categories(self, eggnog_results_df) -> dict:
        eggnog_results_df["COG_category"] = eggnog_results_df["COG_category"].apply(
            lambda x: self.turn_category_to_list(x)
        )
        return eggnog_results_df["COG_category"].to_dict()

    def parse_eggnog_raw_results(self):
        # Load reference proteins
        self.get_reference_fasta_file()
        ref_data = self.load_all_reference_proteins()
        # Load eggnog annotation
        eggnog_results_f = self.out_dir / (self.ref + "_eggNOG_results.csv")
        eggnog_results_df = pd.read_csv(eggnog_results_f, sep="\t", index_col=0, comment="#", names=self.eggnog_headers)
        eggnog_annotation = self.get_assigned_cog_categories(eggnog_results_df)
        # Update reference data 
        ref_data.update(eggnog_annotation)
        self.proteome_cog = ref_data

    def get_protein_subsets(self):
        cols = self.core_protein_table.columns
        proteins = {}
        for col in cols:
            mask = self.core_protein_table[col] == 1
            if "fingerprint" in col:
                proteins["fingerprints"] = self.core_protein_table.loc[mask].index.tolist()
            else:
                proteins[col] = self.core_protein_table.loc[mask].index.tolist()
        self.protein_subsets = proteins

    def calculate_cog_categories_for_proteins(self):
        # Initialize data
        cog_categories = {}
        for category in self.category_annotation:
            cog_categories[category] =  {"proteome":0}
            for protein_subset, proteins in self.protein_subsets.items():
                cog_categories[category][protein_subset] = 0

        # Add proteome data
        for protein, protein_categories in self.proteome_cog.items():
            for category in protein_categories:
                cog_categories[category]["proteome"] += 1

        # Add data for each subset (core, fingerprints)
        for protein_subset, proteins in self.protein_subsets.items():
            for protein in proteins:
                protein_categories = self.proteome_cog[protein] # WIll always exist as key
                for category in protein_categories:
                    cog_categories[category][protein_subset] += 1

        self.cog_categories_data = cog_categories

    def perform_hypergeom_test(self, population_size: int, population_successes: int, sample_size: int, sample_successes: int) -> Tuple:
        expected = (sample_size * population_successes) / population_size
        # p-value
        pval = hypergeom.sf(sample_successes - 1, population_size, population_successes, sample_size)
        if sample_successes < expected or sample_successes == 0:
            pval = 1-hypergeom.sf(sample_successes, population_size, population_successes, sample_size)
        # fold change
        ratio_population = population_successes / population_size
        ratio_sample = sample_successes / sample_size
        if ratio_population == 0 or ratio_sample == 0:
            return pval, 0
        fold_change = ratio_sample / ratio_population
        return pval, fold_change

    def compare_sets_with_hypergeometric_test(self) -> dict:
        background_set = "proteome"
        population_size = len(self.proteome_cog.keys())# Number of proteins in proteome
        dfs = []
        for protein_set, protein_list in self.protein_subsets.items():
            data = {self.ref: {}}
            sample_size = len(protein_list)
            for category in self.cog_categories_data:
                population_successes = self.cog_categories_data[category][background_set]
                sample_successes = self.cog_categories_data[category][protein_set]
                pvalue, fold_change = self.perform_hypergeom_test(population_size, population_successes, sample_size, sample_successes)
                data[self.ref][f"{category}_pvalue"] = pvalue
                data[self.ref][f"{category}_fold_change"] =  fold_change
                # data[self.ref][category] = {"Annotation": self.category_annotation.get(category, "X"), "p-value": pvalue, "Fold_change": fold_change}
            df = pd.DataFrame.from_dict(data, orient="index")
            dfs.append(df)
        return {"core": dfs[0], "fingerprint":dfs[1]}
    
    def highlight_pvalue_on_output(self, row: pd.Series) -> List:
        over_rep_colour = "background-color:green"
        under_rep_colour = "background-color:red"
        return_value = ["" for _ in row.index]
        pvalue_idx = [idx for idx in row.index if "pvalue" in idx]
        fold_change_idx = [idx for idx in row.index if "fold_change" in idx]
        for pvalue_col, fold_change_col in zip(pvalue_idx, fold_change_idx):
            pvalue = row[pvalue_col]
            fold_change = row[fold_change_col]
            pvalue_col_idx = row.index.get_loc(pvalue_col)
            fold_change_col_idx = row.index.get_loc(fold_change_col)
            if type(pvalue) == str:
                return return_value
            if pvalue <= 0.05:
                if fold_change > 1:
                    return_value[pvalue_col_idx] = over_rep_colour
                    return_value[fold_change_col_idx] = over_rep_colour
                if fold_change < 1:
                    return_value[pvalue_col_idx] = under_rep_colour
                    return_value[fold_change_col_idx] = under_rep_colour
        return return_value

    def write_hypergeometric_dfs(self):
        fout = self.out_dir /  "eggNOG_hypergeometric.xlsx"
        excel_writer = pd.ExcelWriter(fout)
        for protein_set, df in zip(self.protein_subsets, self.hypergeometric_dfs):
            df = df.style.apply(self.highlight_pvalue_on_output, axis=1)
            df.to_excel(excel_writer, sheet_name=protein_set)
        excel_writer.close()

    def gather_eggnog_results(self):
        print(f"Parsing results: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        total_data = {"core": [], "fingerprint": []}
        for core_protein_table_f in self.core_protein_table_file_list:
            self.core_protein_table = pd.read_excel(core_protein_table_f, index_col=0)
            self.ref = str(self.core_protein_table.index.name)
            self.parse_eggnog_raw_results()
            self.get_protein_subsets()
            self.calculate_cog_categories_for_proteins()
            tmp_data = self.compare_sets_with_hypergeometric_test()
            total_data["core"].append(tmp_data["core"])
            total_data["fingerprint"].append(tmp_data["fingerprint"])
        self.hypergeometric_dfs = [
            pd.concat(total_data["core"]),
            pd.concat(total_data["fingerprint"])
        ]
        print(f"Writing results: {datetime.now().strftime('%m/%d/%Y, %H:%M:%S')}")
        self.write_hypergeometric_dfs()
    
