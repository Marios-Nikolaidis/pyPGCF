from pathlib import Path
from tqdm import tqdm
import os
from pathlib import Path
import asyncio

class Downloader():
    def __init__(self, taxon: str, assembly_levels: list, out_dir: Path, genes: bool, genomes:bool, protein:bool):
        self.taxon = taxon
        self.assembly_levels = assembly_levels
        self.out_dir = out_dir
        filetypes = []
        if genes:
            filetypes.append("cds")
        if protein:
            filetypes.append("protein")
        if genomes:
            filetypes.append("genome")
        self.filetypes = filetypes

    def download_dehydrated(self):
        cmd = " ".join([
            "datasets download genome taxon",
            self.taxon,
            "--assembly-level",
            ",".join(self.assembly_levels),
            "--include",
            ",".join(self.filetypes),
            f"--assembly-source RefSeq --dehydrated --filename {self.out_dir}/dataset.zip"
            ])
        os.system(cmd)

    def unzip_dehydrated_archive(self):
        cmd = f"unzip {self.out_dir}/dataset.zip -d {self.out_dir}"
        os.system(cmd)

    async def rehydrate_download(self):
        cmd = f"datasets rehydrate --directory {self.out_dir}"
        os.system(cmd)

    def create_output_directories(self):
        outdirs = []
        if "cds" in self.filetypes:
            outdirs.append(self.out_dir/"Genes_fasta_files")
        if "protein" in self.filetypes:
            outdirs.append(self.out_dir/"Protein_fasta_files")
        if "genome" in self.filetypes:
            outdirs.append(self.out_dir/"Genomes_fasta_files")

        for outdir in outdirs:
            outdir.mkdir(exist_ok=True, parents=True)

    def organize_files(self):
        directory = self.out_dir / "ncbi_dataset" / "data"
        genome_id_directories = list(directory.glob("*"))
        for genome_directory in tqdm(genome_id_directories, desc="Organizing files", ascii=True):
            files = genome_directory.glob("*")
            for f in files:
                if "genomic" in f.stem:
                    target = self.out_dir / "Genomes_fasta_files"
                    fout = target / (genome_directory.name + ".fna")
                    f.rename(fout)
                    continue
                if "protein" in f.stem:
                    target = self.out_dir / "Protein_fasta_files"
                    fout = target / (genome_directory.name + ".faa")
                    f.rename(fout)
                    continue
                if "cds_from_genomic" in f.stem:
                    target = self.out_dir / "Genes_fasta_files"
                    fout = target / (genome_directory.name + ".fna")
                    f.rename(fout)
                    continue

    def recursive_directory_search(self, directory: Path):
        items = directory.glob("*")
        for item in items:
            if item.is_dir():
                self.recursive_directory_search(item)
                item.rmdir()
            else:
                item.unlink()

    def remove_dataset_archive(self):
        archive = self.out_dir / "dataset.zip"
        archive.unlink()
        readme = self.out_dir / "README.md"
        readme.unlink()
        ncbi_directory = self.out_dir / "ncbi_dataset"
        self.recursive_directory_search(ncbi_directory)
        ncbi_directory.rmdir()



    async def download(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)
        self.download_dehydrated()
        self.unzip_dehydrated_archive()
        await self.rehydrate_download()
        self.create_output_directories()
        self.organize_files()
        self.remove_dataset_archive()
        
def main():
    taxon = str(810)
    assembly_levels = ["complete", "chromosome"]
    out_dir = Path("/home/marios/Downloads/chlamydia_data")
    genes = False
    genomes = True
    protein = True
    downloader = Downloader(taxon, assembly_levels, out_dir, genes, genomes, protein)
    asyncio.run(downloader.download())

if __name__ == "__main__":
    main()
