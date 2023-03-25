
def download(txid: str):
    cmd = f"datasets download genome taxon {txid} --assembly-level chromosome,complete --assembly-source RefSeq"
    os.system(cmd)
