import requests
from Bio import SeqIO
import json
import asyncio
import time
from pathlib import Path
from tqdm import tqdm
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

class smBGCWeb():
    """
    The base links from
    """
    def __init__(self, out_dir: Path):
        self.api_url = "https://antismash.secondarymetabolites.org/api/v1.0/"
        self.url ="https://antismash.secondarymetabolites.org"
        self.out_dir = out_dir / "smBGC"

    def setup_directories(self):
        self.out_dir.mkdir(exist_ok=True, parents=True)

    def check_status(self, token:str):
        status_url = self.api_url + "status/" + token
        res = requests.get(status_url)
        res_json = res.json()
        state = res_json["state"] # done is finished
        if state == "done":
            return True
        return False

    def keep_track_of_pending(self, tracking_dict: dict, limit: int):
        print(f"I am inside keep track func")
        num_pending = len(tracking_dict)
        while True:
            for key, val in tracking_dict.items():
                print(f"Checking {key}, limit: {limit}")
                if val["done"]: # no need to check. It is already done
                    print("OK")
                    num_pending -= 1
                    continue
                response = val["res"]
                done = self.check_status(response.json()["id"])
                if done:
                    print(f"{key} is now done")
                    tracking_dict[key]["done"] = done
                time.sleep(1)
            if num_pending < limit:
                break

    def download(self, token:str, out_dir:Path):
        status_url = self.api_url + "status/" + token
        print(status_url)
        res = requests.get(status_url)
        res_json = res.json()
        filename = res_json["filename"]
        filename = ".".join(filename.split(".")[:-1]) + ".zip"
        result_url =  self.url + res_json["result_url"].replace("index.html", filename)
        fout = out_dir / filename
        fout.write_bytes(requests.get(result_url).content)

    def post(self, record):
        post_url = self.api_url + "submit"
        payload = {
           "knownclusterblast":"true",
           "clusterblast":"false",
           "subclusterblast":"false",
           "cc_mibig":'false',
           "asf":"false",
           "rre":"false",
           "clusterhmmer":"false",
           "pfam2go":"false",
           "tfbs":"false",
           "jobtype":"antismash7",
           "hmmdetection_strictness":"strict",
           "genefinder":"prodigal",
           "ncbi":record.id
        }
        res = requests.post(post_url, data=payload)
        return res
        
    def load_genome_file(self, genome_file) -> dict:
        parser = SeqIO.parse(str(genome_file), "fasta")
        tracking = {}
        for record in parser:
            res = self.post(record)
            tracking[record.id] = {
                "genome": genome_file.stem,
                "res": res,
                "done": False
            }
        return tracking

if __name__ == "__main__":
    web_mod = smBGCWeb(Path("data/"))
    web_mod.setup_directories()
    input_files = list(Path("data/Genome_fasta_files").glob("*.fna"))
    tracking_dict = {}
    num_pending = 0
    for idx, input_file in enumerate(input_files):
        print(f"Init {idx}, {input_file}")
        file_tracking_dict = web_mod.load_genome_file(input_file)
        num_pending += len(file_tracking_dict.keys())
        tracking_dict.update(file_tracking_dict)
        web_mod.keep_track_of_pending(tracking_dict, 5)
        print(f"Finished {idx}, {input_file}")
    web_mod.keep_track_of_pending(tracking_dict, 1)
    print("Finished everything")
    
    for record, values in tqdm(tracking_dict.items(), desc="Downloading.."):
        genome = values["genome"]
        response = values["res"]
        genome_directory: Path = web_mod.out_dir / genome
        genome_directory.mkdir(exist_ok=True)
        token = response.json()["id"]
        web_mod.download(token, genome_directory)





# class smBGCLocal():
#     def __init__(self, indir: Path, res_dir: Path, cores: int):
#         pass
#         # run_antismash $indir/$file $outdir --cpu 5 --hmmdetection-strictness strict --genefinding-tool prodigal --cb-knownclusters
#
# class smBGCParser():
#     def __init__(self, res_dir: Path, cores: int):
#         self.res_dir = res_dir
#         self.antismash_subdirs = list(self.res_dir.glob("*"))
#         self.cores = cores
#
#     def load_json(self, f: str) -> dict:
#         json_str = ""
#         with open(f, "r") as rf:
#             for line in rf:
#                 line = line.rstrip()
#                 json_str += line
#         return json.loads(json_str)
#
#     def parse_strain_results(self, strain_subdir: Path) -> list:
#         genome = strain_subdir.name
#         json_f  = strain_subdir / (genome + ".json")
#         json_str = self.load_json(str(json_f))
#         data = []
#         for seq_idx, seq_obj in enumerate(json_str["records"]):
#             chrom_id = seq_obj["id"]
#             known_cluster_res = seq_obj["modules"]["antismash.modules.clusterblast"]["knowncluster"]["results"]
#             seq_type = "chromosome"
#             if seq_idx != 0:
#                 seq_type = "plasmid"
#             for idx, region in enumerate(seq_obj["areas"]):
#                 # smbgcs
#                 region_s = region["start"]
#                 region_e = region["end"]
#                 region_p = ";".join(region["products"])
#                 # known cluster blast results
#                 if known_cluster_res is None or len(known_cluster_res[idx]["ranking"]) == 0:
#                     num_hits = 0
#                     description = "X"
#                     database_total_genes = 0
#                     perc_identity = 0
#                 else:
#                     num_hits = known_cluster_res[idx]["ranking"][0][1]["hits"]
#                     database_total_genes = len(known_cluster_res[idx]["ranking"][0][0]["proteins"])
#                     perc_identity = round(num_hits/database_total_genes * 100, 0)
#                     if perc_identity > 100: # There may be a gene duplicate in the query smBGC
#                         perc_identity = 100
#                     description = known_cluster_res[idx]["ranking"][0][0]["description"]
#                 tmp_data = {
#                         "Genome": genome,
#                         "Region": chrom_id,
#                         "Region_s": region_s,
#                         "Region_end": region_e,
#                         "Products": region_p,
#                         "Desc": description,
#                         "Perc identity": perc_identity,
#                         "Source":seq_type
#                 }
#                 data.append(tmp_data)
#         return data
#
#     def write_antismash_results(self):
#         base_name = "antiSMASH_results.xlsx"
#         fout = self.res_dir / base_name
#         self.final_df.to_csv(fout)
#
#     def gather_total_results(self):
#         total_data = []
#         with ProcessPoolExecutor(self.cores) as executor:
#             for data in executor.map(self.parse_strain_results, self.antismash_subdirs):
#                 total_data.append(data)
#         self.final_df = pd.DataFrame(total_data)
#         self.write_antismash_results()
#
