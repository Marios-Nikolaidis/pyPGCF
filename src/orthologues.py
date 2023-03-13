import os
import re
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm
import logging
from pathlib import Path
import numpy as np

class Core_genome_instance():
    def __init__(self, fasta_dir: Path, out_dir: Path, ref: str, 
                 blast_threads: int, blast_evalue: float, blast_method: str,
                 system: str):
        self.ref = ref
        self.out_dir = out_dir
        self.system = system
        self.blast_threads = blast_threads
        self.blast_evalue = blast_evalue
        self.fasta_files = list(fasta_dir.glob("*"))
        self.blast_method = blast_method
        # What was the purpose of the update_fasta?
        # Need to remove the ref from the fasta files
        for x in range(len(self.fasta_files)):
            if self.fasta_files[x].name == self.ref:
                break
        self.fasta_files.pop(x)

    def _get_blast_binaries(self) -> list:
        if self.system == 'linux':
            resources_dir = Path('../resources/Linux/')
            makeblastdb_bin = Path(__file__)/resources_dir/"makeblastdb" #/ .replace(os.path.basename(__file__) + "/", "")
            blastn_bin = Path(__file__)/resources_dir/"blastn" #/ .replace(os.path.basename(__file__) + "/", "")
            diamond_bin = Path(__file__)/resources_dir/"diamond" #/ .replace(os.path.basename(__file__) + "/", "")
        else:
            raise Exception("Only linux platforms are supported at the moment..")

        print(f"""
        makeblastdb_bin: {makeblastdb_bin}
        blastn_bin: {blastn_bin}
        blastp_bin: {blastp_bin}
        diamond_bin: {diamond_bin}
              """)
        if self.blast_method == "blastn":
            blast_bin = blastn_bin
            blast_db_bin = makeblastdb_bin
        elif self.blast_method == "diamond":
            blast_bin = diamond_bin
            blast_db_bin = None
        else:
            raise Exception("Homology search method not supported")
        return [blast_bin, blast_db_bin]

    def _execute_cmd(cmd: str) -> None:
        """Generic function to execute a command"""
        # add logging
        #logging.info(f"Executing command: {cmd}")
        os.system(cmd)

    def _create_blast_db(self, fasta_file: Path, blast_db_bin: Path, 
                         threads: int = 4
                         ) -> database_f: Path :
        """
        Create blast database for a fasta file
        """
        database_f = self.out_dir / "Blast_DB_dir" / fasta_file.name
        if self.blast_method == "diamond":
            cmd = f'{blast_db_bin} makedb --in {fasta_file} --db {database_f} --threads {threads}'
        elif self.blast_method == "blastp":
            cmd = f"{blast_db_bin} -in {fasta_file} -dbtype prot -out {database_f}"
        elif self.blast_method == "blastn":
            cmd = f"{blast_db_bin} -in {fasta_file} -dbtype nucl -out {database_f}"
        _execute_cmd(cmd)
        return database_f

    def _create_blast_cmd(self, blast_bin: Path, fasta_file: Path, database_f: Path, out_file: Path) -> str:
        if self.blast_method == "diamond":
            cmd = f"{blast_bin} blastp --query {fasta_file} --db {database_f} --outfmt 6 --out {out_file} --evalue {self.blast_evalue} --threads {self.blast_threads}"
        else:
            cmd = f"{blast_bin} -query {fasta_file} -db {database_f} -outfmt 6 -out {out_file} -evalue {self.blast_evalue} -num_threads {self.blast_threads}"
        return cmd

    def reciprocal_blast(self):
        """
        """
        # Create the refseq blast database
        blast_bin, blast_db_bin = self._get_blast_binaries()
        ref_fasta = self.in_dir / self.ref
        ref_db = self._create_blast_db(ref_fasta, blast_db_bin, self.out_dir)

        for fasta_file in tqdm(fasta_files, ascii=True, leave=False, desc="Performing reciprocal BLAST"):
            # Create the blast database for the fasta file
            database_f = self._create_blast_db(fasta_file, blastdb_bin, self.out_dir)
            # Perform the blast one way
            fout_f = self.out_dir / "Blastp_output" / (self.ref + "_vs_" + fasta_file.name + ".txt"),
            fout_r = self.out_dir / "Blastp_output" / (fasta_file.name + "_vs_" + self.ref + ".txt"),
            out=out_dir / "Blastp_output" / (self.ref + "_vs_" + fasta_file.name + ".txt"),
            forward_blast = self._create_blast_cmd(ref_fasta, blast_bin, database_f, fout)
            reverse_blast = self._create_blast_cmd(fasta_file, blast_bin, ref_db, fout)
            _execute_cmd(forward_blast)
            _execute_cmd(reverse_blast)

        #    fasta_db = NcbimakeblastdbCommandline(cmd=makeblastdb_bin,
        #                                          dbtype="prot",
        #                                          input_file= str(fasta_file),
        #                                          out=out_dir / "Blast_DB_dir" / fasta_file.name,
        #                                          )
        #    blastp_one_way = NcbiblastnCommandline(cmd=blastp_bin,
        #                                            query= str(fasta_file),
        #                                            db=out_dir / "Blast_DB_dir" / self.ref,
        #                                            outfmt=6,
        #                                            out=out_dir / "Blastp_output" / (fasta_file.name + "_vs_" + self.ref + ".txt"),
        #                                            evalue=evalue,
        #                                            num_threads = threads
        #                                            )
    
        #    blastp_reverse = query=input_dir / "Fasta_files" / self.ref,
        #                                            db=out_dir / "Blast_DB_dir" / fasta_file.name,
        #                                            outfmt=6,
        #                                            out=out_dir / "Blastp_output" / (self.ref + "_vs_" + fasta_file.name + ".txt"),
        #                                            evalue=evalue,
        #                                            num_threads= threads
        #                                            )
        #    # Blast output
        #    # qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore
        #    fasta_db()
        #    blastp_one_way()
        #    blastp_reverse()
        #logging.debug(" Finished reciprocal blast ")
        #return None

    def parse_blast_results(self):
        """
        Reads the txt blast output and create the reciprocal table
        First read the files and filter the redundant information for each protein hit, keep the best hit
        Create the reciprocal blast matrices and write into files
        Input: The ref and fasta files
        Output: Reciprocal files at Working_dir/Reciprocal_files
        Return: None
        """
        def _filterBestEvalues(df):
            """
            Helper function to get the best evalue for each Query protein based on evalue
            Input: Pandas Data frame
            Return: Filtered Pandas data frame
            """
        # Note: Significant speed improvement
            indeces_to_keep = []
            first_col = df.columns[0]
            for grp in df.groupby(first_col):
                indeces_to_keep.append(grp[1].index[0])
            final_df = df.loc[indeces_to_keep]
            final_df = final_df.reset_index()
            final_df = final_df.drop("index",axis=1)
            return final_df
    
        def _pidentFilter(df, SD):
            """
            Helper function to create the distribution of percent identities for all the genes in the data frame
            Filter the hits that are not in the +-2 SD range (user defined)
            """
            pid_std  = df.std( axis=0, numeric_only=True).Pident
            pid_mean = df.mean(axis=0, numeric_only=True).Pident
            drop_condition = pid_mean - (SD * pid_std)
            return df.drop(df[df.Pident < drop_condition].index)
    
        def _reciprocal_Matrix(ref_vs_query_df, query_vs_ref_df):
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
    
        logging.debug(" Parsing BLAST results")
        fasta_files = self.fasta_files
        SD = self.params['SD']
        print("Parsing blast output")
    
        for fasta_file in tqdm(fasta_files, ascii=True):
            if fasta_file.name == self.ref:
                continue
            
            ref_vs_query_file = self.out_dir / "Blastp_output" / (self.ref + "_vs_" + fasta_file.name + ".txt")
            query_vs_ref_file = self.out_dir / "Blastp_output" / (fasta_file.name + "_vs_" + self.ref + ".txt")
    
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
            tmp_ref_vs_query_df = _filterBestEvalues(ref_vs_query_df)
            tmp_query_vs_ref_df = _filterBestEvalues(query_vs_ref_df)
            
            reciprocal_df = _reciprocal_Matrix(tmp_ref_vs_query_df, tmp_query_vs_ref_df)
            reciprocal_df_filt = _pidentFilter(reciprocal_df, SD)

            # Create reciprocal files with headers
            rec_fout = self.out_dir / "Reciprocal_files" / (self.ref + "_vs_" + fasta_file.name + ".reciprocal_file")
            reciprocal_df_filt.to_csv(path_or_buf=rec_fout, sep='\t')
        logging.debug(" Finished parsing results")
    
    def create_cog_matrix(self):
        """
        Create the initial cog matrix with empty (NaN) vectors and replace with the true values
        :return: None
        """
        # Pass all genes of the Refseq into a list
        ref_fasta_f = self.params["in"] / "Fasta_files" / self.ref
        refgenes = SeqIO.parse(ref_fasta_f, "fasta")
        # TODO: Writing the files and then reopening them seems kind of waste of resources
        reciprocal_f_dir = self.out_dir / "Reciprocal_files"
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
        columns = total_cog_matrix.columns

        total_cog_f = self.out_dir / 'totalCOGs.csv'
        total_cog_matrix.to_csv(total_cog_f, sep = '\t', na_rep="X") # Contains all the COGs
        #complete_cog_f = self.out_dir / 'completeCOGs.csv'
        #complete_cog_matrix.to_csv(complete_cog_f, sep = '\t') # Contains the complete COGs
