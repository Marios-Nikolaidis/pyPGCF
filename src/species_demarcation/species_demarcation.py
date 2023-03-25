import os
import numpy as np
import csv
import pathlib
import pandas as pd

def write_tmp_file(files_for_fastani, tmp_file_for_fastani):
    with open(str(tmp_file_for_fastani), "w") as f:
        for file in files_for_fastani:
            f.write(str(file) + "\n")

def perform_fastani(org_list, threads, fraglen, min_fraction, kmer_size, fout):
    os.system("fastANI --ql {} --rl {} -t {} -k {} --fragLen {} minFraction {} -o {}".format(org_list, org_list, threads, kmer_size, fraglen, min_fraction, fout))

def prepare_input_for_mcl(input_file: Path) -> Path:
    """
    This script prepares the input file for fastANI.
    """
    headers = ["query", "target", "ANI", "query_length", "target_length"]
    df = pd.read_csv(input_file, sep="\t", index_col=0, names=headers)
    df["ANI"] = df["ANI"].apply(lambda x: np.nan if x < 95 else x)
    df = df.dropna()
    df = df.drop(columns=["query_length", "target_length"])
    # Write the data to a new file
    fout = input_file.parent / "fastANI_for_mcl.txt"
    df.to_csv(fout, sep="\t", header=False)
    return fout

def run_mcl(fastani_for_mcl, inflation: float, threads: int) -> Path:
    # Make each one a system call
    outdir = fastani_for_mcl.parent
    os.system(f"mcxload -abc {fastani_for_mcl} -o {outdir}/fastANI_mcx_mtrx.txt -write-tab fastANI_annot.tab")
    os.system(f"mcl {outdir}/fastANI_mcx_mtrx.txt -te {threads} -I {inflation} -o {outdir}/fastANI_mcl_out.txt")
    os.system(f"mcxdump -icl {outdir}/fastANI_mcl_out.txt -tabr {outdir}/fastANI_annot.tab -o {outdir}/fastANI_mcx_dump.txt")
    os.system(f"rm {outdir}/fastANI_for_mcl.txt {outdir}/fastANI_mcx_mtrx.txt {outdir}/fastANI_annot.tab {outdir}/fastANI_mcl_out.txt")
    os.system(f"mv {outdir}/fastANI_mcx_dump.txt {outdir}/fastANI_clusters.tsv")
    return outdir/"fastANI_clusters.tsv"


def parse_mcx_output(fastani_from_mcl):
    outdir = fastani_from_mcl.parent
    results = {}
    clust_num = 0
    for lines in csv.reader(open(str(fastni_from_mcl), "r"), delimiter="\t"):
        for l in lines:
            results[l] = clust_num
        clust_num += 1
    df = pd.DataFrame.from_dict(results, orient="index")
    df.columns = ["ClustNum"]
    df["Clust"] = df["ClustNum"].apply(lambda x: "C" + str(x))
    df = df.drop("ClustNum", axis=1)
    fastani_from_mcl.unlink()
    fout = outdir / "fastANI_mcl_clusters.xlsx"
    df.to_excel(fout)

def main(indir: Path, outdir: Path, threads:int, fraglen: int, min_fraction: float, kmer_size: int, mcl_inflation: float):
    # Check if input is file or directory
    files_for_fastani = indir.glob("*")
    tmp_file_for_fastani = outdir / "fastANI_input.txt"
    write_tmp_file(files_for_fastani, tmp_file_for_fastani)
    fastani_out = outdir / "fastANIoutput.txt"
    perform_fastani(tmp_file_for_fastani, threads, fraglen, min_fraction, kmer_size, fastani_out)
    fastani_for_mcl = prepare_input_for_mcl(fastani_out)
    fastani_from_mcl = run_mcl(fastani_for_mcl, mcl_inflation, threads)
    parse_mcx_output(fastani_from_mcl)
