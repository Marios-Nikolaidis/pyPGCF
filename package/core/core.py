""" Calculate the core proteins and fingerprints based on a given reference genome"""
import pandas as pd
from typing import List
from pathlib import Path

def split_genomes_into_groups(df: pd.DataFrame, ref: str, species_df,  omplete:bool = False) -> List:
    """
    Split the genomes into groups based on the species cluster (cluster_col)
    If complete is True, then the complete table (taxon) is used
    """
    if complete:
        group_orgs = list(df.index)
        group_orgs.remove(ref)
        return [group_orgs, []]
    ref_cluster = species_df.loc[ref].values[0]
    group_orgs = list(df[df[cluster_col] == ref_cluster].index)
    non_group_orgs = list(df.index)
    for org in group_orgs:
        non_group_orgs.remove(org)
    if ref in group_orgs:
        group_orgs.remove(ref)
    return [group_orgs, non_group_orgs]


def calculate_protein_presence(
    df: pd.DataFrame, org_list: List[List[str]]
) -> pd.DataFrame:
    group_orgs, non_group_orgs = org_list
    tmpdf = df.copy().drop(df.columns, axis=1)
    tmpdf["Group orthologues"] = df[group_orgs].count(axis=1)
    tmpdf["Group orthologues"] = tmpdf["Group orthologues"] + 1
    tmpdf["Group orthologues %"] = round(
        (tmpdf["Group orthologues"] / (len(group_orgs) + 1)) * 100, 2
    )
    tmpdf["Non group orthologues"] = df[non_group_orgs].count(axis=1)
    tmpdf["Non group orthologues %"] = round(
        tmpdf["Non group orthologues"] / len(non_group_orgs) * 100, 2
    )
    return tmpdf

def identify_core_proteins(df: pd.DataFrame, core_perc: float) -> List:
    """"""
    core_proteins = df[df["Group orthologues %"] >= core_perc].index.to_list()
    fingerprints = df[df[(df["Group orthologues %"] == 100) & (df["Non group orthologues %"] == 0)]].index.to_list()
    data = {protein: {f"Core_{core_perc}%": 1, "Fingerprint": 0} for protein in core_proteins}
    for protein in fingerprints:
        data[p]["Fingerprint"] = 1
    return data

def main(ref: str, cluster_col: str, cogs_fin: Path, species_file: Path, core_perc: float, outdir: Path) -> str:
    """
    """
    outdir.mkdir(exist_ok=True)
    # Read the COGs matrix
    df = pd.read_csv(cogs_fin, sep="\t", index_col=0)
    # Split the genomes into groups
    complete = True
    species_df = None
    if species_file != None:
        species_df = pd.read_excel(species_file, index_col=0)
        complete = False
    org_list = split_genomes_into_groups(df, ref, species_df, complete=complete)
    # Calculate the presence of each protein
    final_df = calculate_protein_presence(df, org_list)
    # Identify the core proteins
    core_protein_data = identify_core_proteins(final_df, core_perc)
    core_protein_data_df = pd.DataFrame.from_dict(core_protein_data, orient="index")
    if complete:
        fout = outdir / f"{ref}_core_proteins_{core_perc}.csv"
    else:
        fout = outdir / f"{ref}_core_proteins_{core_perc}_species.csv"
    core_protein_data_df.to_csv(outdir / f"{ref}_core_proteins_{core_perc}.csv")
    return fout
