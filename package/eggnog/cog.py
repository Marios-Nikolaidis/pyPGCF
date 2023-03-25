"""Add functional annotation to the proteins in the core genes excels"""
import pandas as pd
from collections import Counter
from Bio import SeqIO
import pathlib

def load_all_reference_proteins(fasta_file: Path) -> dict:
    record = SeqIO.parse(str(fasta_file), "fasta")
    return {protein.id: ["S"] for protein in record}

def get_category(x: str) -> str:
    if x == "-":
        return ["S"]
    return list(x)

def load_eggnog_results(eggnog_results: Path) -> dict:
    headers = [ "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", 
        "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
        "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC",
        "CAZy", "BiGG_Reaction", "PFAMs"
    ]
    annotation_f = "local_eggnog_output/rep_proteomes_eggnog.emapper.annotations"
    annotation_df = pd.read_csv(eggnog_results, sep="\t", index_col=0, comment="#", names=headers)
    annotation_df["COG_category"] = annotation_df["COG_category"].apply(
        lambda x: get_category(x)
    )
    return annotation_df["COG_category"].to_dict()

def load_protein_list(protein_list: Path) -> list:
    """
    """
    data = pd.read_excel(protein_list, index_col=0)
    num_cols = df.shape[1]
    proteins = [{}, {}]
    if num_cols == 1: # Only core proteins
        col = data.columns[0]
        proteins[0] = {p: "S" for prot in data[data[col] == 1].index}
    elif num_cols == 2: # Core proteins and fingerprints
        col = data.columns[1]
        proteins[1] = {p: "S" for prot in data[data[col] == 1].index}
    return proteins


def main(fasta_file_in: Path, eggnog_results: Path, protein_list: Path, output_dir: Path):
    category_annotation = {
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
    # Load reference proteins
    ref_data = load_all_reference_proteins(fasta_file_in)
    # Load eggnog annotation
    annotation_df = load_eggnog_results(eggnog_results)
    # Update reference data 
    ref_data.update(annotation_df)
    # Load protein list
    genes_data = load_protein_list(protein_list)
     
    
    #fout = "Genome_functional_categories_20221117.xlsx"
    #data_df = pd.DataFrame.from_dict(data_dictionary, orient="index")
    #data_df.to_excel(fout)
