import os
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import json
import pandas as pd
from tqdm import tqdm
#import requests
#from Bio import SeqIO
#import asyncio
#import time
#from tqdm import tqdm

class smBGCInstaller():
    def __init__(self, debug: bool=False):
        self.debug = debug

    def install_antismash(self):
        print("Installing antiSMASH... this will take some time")
        cmds = [
                #"conda config --add channels defaults",
                #"conda config --add channels bioconda",
                #"conda config --add channels conda-forge",
            #"conda config --set channel_priority strict",
            #"conda install -y -c bioconda antismash"
            "conda install -y antismash"
        ]
        for cmd in cmds:
            os.system(cmd)

    def install_databases(self):
        print("Downloading databases of antiSMASH")
        cmd = "download-antismash-databases"
        os.system(cmd)

if __name__ == "__main__":
    installer = smBGCInstaller()
    installer.install_antismash()

class smBGCLocalRunner():
    def __init__(self, genome_fasta_dir: Path, out_dir: Path, cores: int, strictness: str, genefinding_tool: str):
        self.genome_fasta_dir = genome_fasta_dir
        self.antismash_raw_results_dir = out_dir / "smBGC" / "antismash_raw_results"
        self.cores = cores
        self.strictness = strictness
        self.genefinding_tool = genefinding_tool

    def create_antismash_cmd(self, genome_fasta_file: Path, genome_out_dir: Path):
        cmd = " ".join([
            "antismash",
            f"--output-dir {genome_out_dir}",
            f"--cpu {self.cores}",
            f"--hmmdetection-strictness {self.strictness}",
            f"--genefinding-tool {self.genefinding_tool}",
            "--cb-knownclusters",
            f"{genome_fasta_file.absolute()}",
            ])
        return cmd

    def run_antismash_cmd(self, cmd: str):
        os.system(cmd)

    def analyze_genomes(self):
        genome_files = list(self.genome_fasta_dir.glob("*"))
        for genome_file in tqdm(genome_files, desc="Running antiSMASH", ascii=True, leave=True):
            print(genome_file)
            genome_out_dir = self.antismash_raw_results_dir / genome_file.stem
            genome_out_dir.mkdir(exist_ok=True, parents=True)
            cmd = self.create_antismash_cmd(genome_file, genome_out_dir)
            self.run_antismash_cmd(cmd)

class smBGCParser():
    def __init__(self, out_dir: Path, cores: int):
        self.out_dir = out_dir / "smBGC"
        self.antismash_raw_result_dir = self.out_dir / "antismash_raw_results"
        self.antismash_subdirs = list(self.antismash_raw_result_dir.glob("*"))
        self.parsing_cores = cores

    def load_json(self, f: str) -> dict:
        json_str = ""
        with open(f, "r") as rf:
            for line in rf:
                line = line.rstrip()
                json_str += line
        return json.loads(json_str)

    def get_known_cluster_results(self, seq_obj):
        clusterblast_res = seq_obj["modules"]["antismash.modules.clusterblast"]["knowncluster"]["results"]
        return clusterblast_res

    def has_identified_smbgcs(self, identified_smbgcs: list) -> bool:
        if len(identified_smbgcs) == 0:
            return False
        return True

    def parse_strain_results(self, strain_subdir: Path) -> pd.DataFrame:
        genome = strain_subdir.name
        json_f  = strain_subdir / (genome + ".json")
        json_str = self.load_json(str(json_f))
        data = []
        for seq_idx, seq_obj in enumerate(json_str["records"]):
            identified_smbgcs = seq_obj["areas"]
            has_identified_smbgcs = self.has_identified_smbgcs(identified_smbgcs)
            if not has_identified_smbgcs:
                continue
            known_cluster_results = self.get_known_cluster_results(seq_obj)
            seq_type = "chromosome"
            if seq_idx != 0:
                seq_type = "plasmid"
            region_id = seq_obj["id"]
            for idx, area in enumerate(identified_smbgcs):
                # smbgcs
                region_s = area["start"]
                region_e = area["end"]
                region_p = ";".join(area["products"])
                # known cluster blast results
                if known_cluster_results is None or len(known_cluster_results[idx]["ranking"]) == 0:
                    num_hits = 0
                    description = "X"
                    database_total_genes = 0
                    perc_identity = 0
                else:
                    num_hits = known_cluster_results[idx]["ranking"][0][1]["hits"]
                    database_total_genes = len(known_cluster_results[idx]["ranking"][0][0]["proteins"])
                    perc_identity = round(num_hits/database_total_genes * 100, 0)
                    if perc_identity > 100: # There may be a gene duplicate in the query smBGC
                        perc_identity = 100
                    description = known_cluster_results[idx]["ranking"][0][0]["description"]
                tmp_data = {
                        "Genome": genome,
                        "Source":seq_type,
                        "Region": region_id,
                        "Region_s": region_s,
                        "Region_end": region_e,
                        "Products": region_p,
                        "Most similar known cluster": description,
                        "Perc. similarity": perc_identity,
                }
                data.append(tmp_data)
        return pd.DataFrame(data)

    def write_antismash_results(self):
        excel_filename = "antiSMASH_results.xlsx"
        fout = self.out_dir / excel_filename
        self.final_df.to_excel(fout, index=False)

    def gather_results(self):
        total_data = []
        with ProcessPoolExecutor(self.parsing_cores) as executor:
            for dataframe in tqdm(executor.map(self.parse_strain_results, self.antismash_subdirs), total=len(self.antismash_subdirs), desc="Parsing antiSMASH results", ascii=True, leave=True):
                total_data.append(dataframe)
        self.final_df = pd.concat(total_data)
        self.write_antismash_results()

