import os
import numpy as np
import logging
import csv
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from typing import Union, Generator, List

class SpeciesDemarcator():
    def __init__(self,  in_dir: Path, out_dir: Path, fastani_cores:int, kmer:int, fraglen:int, minfrac:float, inflation:int, mcl_cores:int):
        self.in_dir = in_dir
        self.out_dir = out_dir / "Species_demarcation"
        self.fastani_cores = fastani_cores
        self.kmer = kmer
        self.fraglen = fraglen
        self.minfrac = minfrac
        self.inflation = inflation
        self.mcl_cores = mcl_cores

    def create_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)

    def create_input_for_fastani(self, files_for_fastani: Union[List, Generator], tmp_file_for_fastani: Path) -> None:
        print("Preparing FastANI input")
        with open(str(tmp_file_for_fastani), "w") as f:
            for file in files_for_fastani:
                f.write(str(file) + "\n")
    
    def perform_fastani(self, org_list: Path, fout: Path) -> None:
        print("Performing FastANI...")
        cmd = "fastANI --ql {} --rl {} -t {} -k {} --fragLen {} --minFraction {} -o {} > /dev/null 2>&1".format(org_list, org_list, self.fastani_cores, self.kmer, self.fraglen, self.minfrac, fout)
        os.system(cmd)
    
    def prepare_input_for_mcl(self, input_file: Path) -> Path:
        """
        This script prepares the input file for fastANI.
        """
        headers = ["query", "target", "ANI", "query_length", "target_length"]
        df = pd.read_csv(input_file, sep="\t", index_col=0, names=headers)
        df["ANI"] = df["ANI"].apply(lambda x: np.nan if x < 95 else x)
        df = df.dropna()
        df = df.drop(columns=["query_length", "target_length"])
        fout = input_file.parent / "fastANI_for_mcl.txt"
        df.to_csv(fout, sep="\t", header=False)
        return fout
    
    def run_mcl(self, fastani_for_mcl: Path) -> Path:
        # Make each one a system call
        print("Running MCL clustering")
        outdir = fastani_for_mcl.parent
        os.system(f"mcxload -abc {fastani_for_mcl} -o {outdir}/fastANI_mcx_mtrx.txt -write-tab {outdir}/fastANI_annot.tab > /dev/null 2>&1")
        os.system(f"mcl {outdir}/fastANI_mcx_mtrx.txt -te {self.mcl_cores} -I {self.inflation} -o {outdir}/fastANI_mcl_out.txt > /dev/null 2>&1")
        os.system(f"mcxdump -icl {outdir}/fastANI_mcl_out.txt -tabr {outdir}/fastANI_annot.tab -o {outdir}/fastANI_mcx_dump.txt > /dev/null 2>&1")
        os.system(f"rm {outdir}/fastANI_for_mcl.txt {outdir}/fastANI_mcx_mtrx.txt {outdir}/fastANI_annot.tab {outdir}/fastANI_mcl_out.txt {outdir}/FastANI_input.txt")
        os.system(f"mv {outdir}/fastANI_mcx_dump.txt {outdir}/fastANI_clusters.tsv")
        print("Done")
        return outdir/"fastANI_clusters.tsv"
    
    def parse_mcx_output(self, fastani_from_mcl: Path) -> None:
        print("Parsing MCL output")
        outdir = fastani_from_mcl.parent
        results = {}
        clust_num = 0
        for lines in csv.reader(open(str(fastani_from_mcl), "r"), delimiter="\t"):
            for l in lines:
                results[l] = clust_num
            clust_num += 1
        df = pd.DataFrame.from_dict(results, orient="index")
        df.columns = ["ClustNum"]
        df["FastANI_species"] = df["ClustNum"].apply(lambda x: "C" + str(x))
        df = df.drop("ClustNum", axis=1)
        df.index = [idx.split("/")[-1] for idx in df.index]
        df.index = [".".join(idx.split(".")[:-1]) for idx in df.index]
        fastani_from_mcl.unlink()
        fout = outdir / "FastANI_species_clusters.xlsx"
        df.to_excel(fout)
    
    def assign_species(self):
        # Check if input is file or directory
        self.create_directories()
        files_for_fastani = self.in_dir.glob("*")
        tmp_file_for_fastani = self.out_dir / "FastANI_input.txt"
        self.create_input_for_fastani(files_for_fastani, tmp_file_for_fastani)
        fastani_out = self.out_dir / "FastANI.tsv"
        self.perform_fastani(tmp_file_for_fastani, fastani_out)
        fastani_for_mcl = self.prepare_input_for_mcl(fastani_out)
        fastani_from_mcl = self.run_mcl(fastani_for_mcl)
        self.parse_mcx_output(fastani_from_mcl)
        print("Done")
