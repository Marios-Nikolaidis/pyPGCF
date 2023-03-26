import os
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import logging
from pathlib import Path
import numpy as np
from typing import Tuple

class Orthologues_identifier():
    def __init__(self, fasta_dir: Path, out_dir: Path, ref: str, 
                 blast_threads: int, blast_evalue: float, blast_method: str,
                 filter: bool):
        self.ref = ref
        self.in_dir = fasta_dir.parent
        self.out_dir = out_dir / "Orthologues"
        self.blast_threads = blast_threads
        self.blast_evalue = blast_evalue
        self.fasta_files = list(fasta_dir.glob("*"))
        self.blast_method = blast_method
        self.threads = blast_threads
        self.filter = filter
        # Remove the ref from the fasta files
        for index, item in enumerate(self.fasta_files):
            if self.ref in item.name:
                self.ref_fasta = item
                self.fasta_files.pop(index)
                break

    def _set_directories(self) -> None:
        directories = [
                self.out_dir / self.ref / "Blast_results",
                self.out_dir / self.ref / "Best_reciprocal_hits",
                self.out_dir / self.ref / "Blast_DB",
                ]
        for directory in directories:
            directory.mkdir(exist_ok=True, parents=True)
        

    def _get_blast_binaries(self) -> None:
        resources_dir = Path('../../resources/')
        makeblastdb_bin = Path(__file__).parent/resources_dir/"makeblastdb" #/ .replace(os.path.basename(__file__) + "/", "")
        blastn_bin = Path(__file__).parent/resources_dir/"blastn" #/ .replace(os.path.basename(__file__) + "/", "")
        diamond_bin = Path(__file__).parent/resources_dir/"diamond" #/ .replace(os.path.basename(__file__) + "/", "")

        if self.blast_method not in ("blastn", "diamond"):
            raise Exception("Homology search method not supported")
        self.blast_bin = diamond_bin
        if self.blast_method == "blastn":
            self.blast_bin = blastn_bin
            self.blast_db_bin = makeblastdb_bin

    def _execute_cmd(self, cmd: str) -> int:
        """Generic function to execute a command"""
        # add logging
        #logging.info(f"Executing command: {cmd}")
        return os.system(cmd)

    def _create_blast_db(self, fasta_file: Path) -> Tuple[int, Path] :
        """
        Create blast database for a fasta file
        """
        database_f = self.out_dir / self.ref / "Blast_DB" / fasta_file.name
        cmd = f"{self.blast_bin} makedb --in {fasta_file} --db {database_f} --threads {self.threads}" # DIAMOND
        if self.blast_method == "blastn":
            cmd = f"{self.blast_db_bin} -in {fasta_file} -dbtype nucl -out {database_f}"
        database_f = database_f.with_suffix(".faa.dmnd")
        return self._execute_cmd(cmd), database_f

    def _create_blast_cmd(self, fasta_file: Path, database_f: Path, out_file: Path) -> str:
        # DIAMOND case first
        cmd = f"{self.blast_bin} blastp --query {fasta_file} --db {database_f} --outfmt 6 --out {out_file} --evalue {self.blast_evalue} --threads {self.blast_threads}"
        if self.blast_method == "blastn":
            cmd = f"{self.blast_bin} -query {fasta_file} -db {database_f} -outfmt 6 -out {out_file} -evalue {self.blast_evalue} -num_threads {self.blast_threads}"
        return cmd

    def setup(self):
        self._set_directories()
        self._get_blast_binaries()

    def reciprocal_blast(self):
        """ """
        # Create the refseq blast database
        _, ref_db = self._create_blast_db(self.ref_fasta)

        for fasta_file in tqdm(self.fasta_files, ascii=True, leave=False, desc="Performing reciprocal BLAST"):
            # Create the blast database for the fasta file
            _, database_f = self._create_blast_db(fasta_file)
            # Perform the blast one way
            fout_f = self.out_dir / self.ref / "Blast_results" / (fasta_file.name + "_vs_reference.txt")
            fout_r = self.out_dir / self.ref / "Blast_results" / (fasta_file.name + "_vs_reference_reverse.txt")
            forward_blast = self._create_blast_cmd(self.ref_fasta, database_f, fout_f)
            reverse_blast = self._create_blast_cmd(fasta_file, ref_db, fout_r)
            self._execute_cmd(forward_blast)
            self._execute_cmd(reverse_blast)

    def parse_blast_results(self):
        """
        Reads the txt blast output and create the reciprocal table
        First read the files and filter the redundant information for each protein hit, keep the best hit
        Create the reciprocal blast matrices and write into files
        Input: The ref and fasta files
        Output: Best reciprocal hits per strain
        Return: None
        """
        def _get_best_subject(df: pd.DataFrame) -> pd.DataFrame:
            """
            Helper function to keep the best subject hit per query
            Input: Pandas Data frame
            Return: Filtered Pandas data frame
            """
            indeces_to_keep = []
            first_col = df.columns[0]
            for grp in df.groupby(first_col):
                indeces_to_keep.append(grp[1].index[0])
            return df.loc[indeces_to_keep].reset_index(drop=True)
    
        def _orthologue_filter(df: pd.DataFrame) -> pd.DataFrame:
            """
            Helper function to create the distribution of percent identities for all the genes in the data frame
            Filter the hits that have perc. identity lower than -2 SD
            """
            pid_std  = df.Pident.std( axis=0, numeric_only=True)
            pid_mean = df.Pident.mean(axis=0, numeric_only=True)
            drop_condition = pid_mean - (2 * pid_std)
            return df.drop(df[df.Pident < drop_condition].index)
    
        def _create_reciprocal_matrix(ref_vs_query_df, query_vs_ref_df):
            """
            Helper function which compare the two dataframes given line by line 
            to identify reciprocal hits
            Input: The one way and the reverse way dataframe
            Return a new data frame which contains only the reciprocal hits
            """
            tmp_df = pd.merge(ref_vs_query_df, query_vs_ref_df, on="RefSeq")
            df = tmp_df.drop(tmp_df[tmp_df.QuerySeq_x != tmp_df.QuerySeq_y].index)
            case1 = df[df['Pident_x'] >= df['Pident_y']]
            case1 = case1[['RefSeq', 'QuerySeq_x', 'Pident_x', 'Evalue_x']]
            case1 = case1.rename(columns={'QuerySeq_x':'QuerySeq', 'Pident_x':'Pident', 'Evalue_x':'Evalue'})
            case2 = df[df['Pident_x'] < df['Pident_y']]
            case2 = case2.rename(columns={'QuerySeq_y':'QuerySeq', 'Pident_y':'Pident', 'Evalue_y':'Evalue'})
            df = pd.concat([case1,case2])
            return df[['RefSeq', 'QuerySeq', 'Pident', 'Evalue']]
    
        logging.debug("Parsing BLAST results")
        fasta_files = self.fasta_files
        print("Parsing blast output")
    
        for fasta_file in tqdm(fasta_files, ascii=True):
            if fasta_file.name == self.ref:
                continue
            ref_vs_query_file = self.out_dir / self.ref / "Blast_results" / (fasta_file.name + "_vs_reference.txt")
            query_vs_ref_file = self.out_dir / self.ref / "Blast_results" / (fasta_file.name + "_vs_reference_reverse.txt")
    
            c1 = ["RefSeq", "QuerySeq", "Pident", "Evalue"]
            c2 = ["QuerySeq", "RefSeq", "Pident", "Evalue"]
            
            ref_vs_query_df = pd.read_csv(ref_vs_query_file,  sep="\t", names=c1, usecols=[0, 1, 2, 10] )
            query_vs_ref_df = pd.read_csv(query_vs_ref_file,  sep="\t", names=c2, usecols=[0, 1, 2, 10] )
            
            ref_vs_query_df.iloc[:, 0] = ref_vs_query_df.iloc[:, 0].astype(str) + "_" + self.ref
            ref_vs_query_df.iloc[:, 1] = ref_vs_query_df.iloc[:, 1].astype(str) + "_" + fasta_file.name
    
            query_vs_ref_df.iloc[:, 0] = query_vs_ref_df.iloc[:, 0].astype(str) + "_" + fasta_file.name
            query_vs_ref_df.iloc[:, 1] = query_vs_ref_df.iloc[:, 1].astype(str) + "_" + self.ref
    
            # Blast output
            # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
            ref_vs_query_df = _get_best_subject(ref_vs_query_df)
            query_vs_ref_df = _get_best_subject(query_vs_ref_df)
            
            reciprocal_df = _create_reciprocal_matrix(ref_vs_query_df, query_vs_ref_df)
            if self.filter:
                reciprocal_df = _orthologue_filter(reciprocal_df)
            # Create reciprocal files with headers
            rec_fout = self.out_dir / self.ref / "Best_reciprocal_hits" / (fasta_file.stem + "_BRH.txt")
            reciprocal_df.to_csv(rec_fout, sep='\t')
        logging.debug(" Finished parsing results")
    
    def create_cog_matrix(self):
        """
        Create the initial cog matrix with empty (NaN) vectors and replace with the true values
        :return: None
        """
        # Pass all genes of the Refseq into a list
        ref_fasta_f = self.in_dir / "Protein_fasta_files" / self.ref
        refgenes = SeqIO.parse(ref_fasta_f, "fasta")
        # TODO: Writing the files and then reopening them seems kind of waste of resources
        reciprocal_f_dir = self.out_dir / self.ref / "Best_reciprocal_hits"
        reciprocal_files = list(reciprocal_f_dir.glob("*"))
    
        def _initCogMatrix(refgenes, rec_files):
            cog_matrix_dict = {gene.id: {} for gene in refgenes}
            for rec_file in rec_files:
                _, query_org_name = rec_file.stem.split("_vs_")
                for gene in cog_matrix_dict:
                    cog_matrix_dict[gene][query_org_name] = np.nan
            return cog_matrix_dict
    
        def _expandCogMatrix(init_cog_matrix_dict, reciprocal_files):
            """
            Replace the NaN values in the Pandas dataframe with the corresponding values
            :return: expanded cogmatrix as pandas dataframe
            """
            for reciprocal_file in reciprocal_files:
                _, query_org_name = reciprocal_file.stem.split("_vs_")
                header = True
                for lines in open(reciprocal_file, mode="r"):
                    if header:
                        header=False
                        continue
                    refseq, queryseq, *_ = lines.rstrip().split("\t")
                    #refseq = refseq.replace('_' + ref,'') Not necessary anymore
                    init_cog_matrix_dict[refseq][query_org_name] = queryseq
            init_cog_matrix = pd.DataFrame.from_dict(init_cog_matrix_dict, orient="index")
            init_cog_matrix.index.name = "RefSeq"
            return init_cog_matrix
    
        init_cog_matrix = _initCogMatrix(refgenes, reciprocal_files)
        total_cog_matrix = _expandCogMatrix(init_cog_matrix, reciprocal_files)
        total_cog_f = self.out_dir / self.ref / "totalCOGs.csv"
        total_cog_matrix.to_csv(total_cog_f, sep = '\t', na_rep="X") # Contains all the COGs
