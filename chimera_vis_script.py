import numpy as np
import ast
import pandas as pd
import os
import re
import itertools
from Bio import Entrez
from chimerax.core.commands import run
from IPython import get_ipython
from IPython.display import display, Image
from pathlib import Path
import requests

# =========================
# Paths
# =========================
folder = os.path.dirname(os.path.abspath(__file__))
print(folder)

png_folder = os.path.join(folder, "2d_map", "png")
pdb_folder = os.path.join(folder, "2d_map", "pdb")
df_file = os.path.join(folder, "output_files", "hits_with_uniprot_rpkm_aco.csv")

Entrez.email = "pde.caltech@gmail.com"


# =========================
# PNG width
# =========================
set_width_png = input("Enter width for PNG display (default 750): ")
set_width_png = int(set_width_png) if set_width_png.isdigit() else 750


# =========================
# Utilities
# =========================
def clear_shell():
    ipy = get_ipython()
    if ipy:
        ipy.system("clear")




def get_residue_ranges(df, accession):
    subset = df[df["Accession"] == accession]
    ranges = []

    for r in subset["Residue Range"].dropna():
        r = str(r).strip()
        if not r or r.lower() in {"nan", "none"}:
            continue

        r = r.replace("–", "-").replace("—", "-")

        # split architectural blocks
        blocks = re.split(r"[;]", r)

        for block in blocks:
            block = block.replace(" ", "")
            matches = re.findall(r"(\d+)-(\d+)", block)
            for start, end in matches:
                start, end = int(start), int(end)
                if start < end:
                    ranges.append((start, end))

    return ranges



def expand_ranges(ranges):
    """Convert [(a,b), (c,d)] → set of all residues."""
    return set(itertools.chain.from_iterable(range(a, b + 1) for a, b in ranges))


def get_alphafold_metadata(accession):
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession}"
    r = requests.get(url, timeout=10)
    r.raise_for_status()
    data = r.json()

    if not data:
        raise ValueError("No AlphaFold data found")

    entry = data[0]
    return {
        "gene_name": entry.get("gene"),
        "organism": entry.get("organismScientificName"),
        "uniprot": entry.get("uniprotAccession")
    }


def search_gene_id(gene_name, organism):
    term = f"{gene_name}[Gene Name] AND {organism}[Organism]"
    handle = Entrez.esearch(db="gene", term=term, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"][0] if record["IdList"] else None







def print_info(accession, extra_info=None):
    try:
        # 1️⃣ Ask AlphaFold FIRST
        af = get_alphafold_metadata(accession)
        gene_name = af["gene_name"]
        organism = af["organism"]

        if not gene_name:
            raise ValueError("AlphaFold did not return a gene name")

        # 2️⃣ Search NCBI Gene properly
        gene_id = search_gene_id(gene_name, organism)
        if not gene_id:
            session.logger.error("No Gene ID found from AlphaFold gene name")
            return

        summary_handle = Entrez.esummary(db="gene", id=gene_id)
        summary_record = Entrez.read(summary_handle)
        summary_handle.close()

        gene_data = summary_record['DocumentSummarySet']['DocumentSummary'][0]

        clear_shell()
        if extra_info:
            session.logger.info(extra_info)

        run(session, "log clear")

        session.logger.info(f"Accession: {accession}")
        session.logger.info(f"Gene ID:   {gene_id}")
        session.logger.info(f"Gene Name: {gene_data['Name']}")
        session.logger.info(f"Description: {gene_data['Description']}")
        session.logger.info("\n--- Function Summary ---")
        session.logger.info(gene_data['Summary'])

        alphafold_link = f"https://alphafold.ebi.ac.uk/entry/{accession}"
        run(session, f"log html \"<a href='{alphafold_link}'>AlphaFold Entry</a>\"")

    except Exception as e:
        session.logger.error(f"An error occurred: {e}")





def gen_png(pdb_folder, df1):
    accession = df1["Accession"].iloc[0]

    files = os.listdir(pdb_folder)
    matches = [f for f in files if f.startswith(f"AF-{accession}")]

    if not matches:
        raise FileNotFoundError(f"No PDB for {accession}")

    pdb_path = os.path.join(pdb_folder, matches[0])


    run(session, "close all")

    run(session, f'open "{pdb_path}"')
    run(session, "set bgColor white")
    run(session, "select all")
    run(session, "color sel #8C7A6f")

    run(session, "select /A:1-3")
    run(session, "color sel magenta")
    run(session, "select clear")

    # --- Serp (RED)
    serp_raw = df1["Serp_Interval"].iloc[0]
    serp_ranges = ast.literal_eval(serp_raw) if isinstance(serp_raw, str) else []
    serp_set = expand_ranges(serp_ranges)

    # --- Domains (BLUE)
    domain_ranges = get_residue_ranges(df1, accession)
    domain_set = expand_ranges(domain_ranges)

    # --- Overlap (GREEN)
    overlap = serp_set & domain_set
    serp_only = serp_set - overlap
    domain_only = domain_set - overlap

    def color_residues(res_set, color):
        if not res_set:
            return
        res = sorted(res_set)
        start = res[0]
        for i in range(1, len(res)):
            if res[i] != res[i - 1] + 1:
                run(session, f"select :{start}-{res[i - 1]}")
                run(session, f"color sel {color}")
                start = res[i]
        run(session, f"select :{start}-{res[-1]}")
        run(session, f"color sel {color}")
        run(session, "select clear")

    color_residues(serp_only, "red")
    color_residues(domain_only, "blue")
    color_residues(overlap, "green")


def render_accession(acc):
    df = pd.read_csv(df_file, dtype=str)
    df = df[df["Accession"] == acc]

    if df.empty:
        raise ValueError(f"{acc} not found in CSV")

    gen_png(pdb_folder, df)


# =========================
# Main loop
# =========================
mode = input("Do you want to load a CSV file with accession numbers? (y/n) (default: n): ").lower()

if mode == "y":

    file = input("Enter path to CSV file with accession numbers: ").strip()

    # remove surrounding single or double quotes if present
    if (file.startswith("'") and file.endswith("'")) or \
    (file.startswith('"') and file.endswith('"')):
        file = file[1:-1]

    p = Path(file)
    stem = p.stem
    start = int(stem[0])
    print('wait')
    if start % 2 == 0:
        start = 'NAC binds after Domain starts'
    else:
        start = 'NAC binds before Domain starts'
    group = stem.split('.', 1)[1]
    group =  group + ' ---> ' + start
    pd_df = pd.read_csv(file)
    accession_column = pd_df['Accession']
    accession_column = pd_df['Accession'].apply(lambda x: x[3:3+6] if len(x) > 6 else x)

    kk = 0
    while True:
        accession = accession_column[kk]
        
        accession = accession.strip().upper()
        png = os.path.join(png_folder, f"{accession}.png")
        '''pdb = os.path.join(pdb_folder, f"{accession}.glb")

        if not os.path.exists(pdb):
            print(f"PDB file for accession {accession} not found at {pdb}.")
            continue'''
        render_accession(accession)
        clear_shell()


        print_info(accession, group)
        display(Image(filename=png, width=set_width_png))
        kk += 1
        if kk >= len(accession_column):
            print("Reached end of accession list.")
            break
        aa = input("Press Enter to go to next accession, or type 'stop' to end: ")

        if aa.lower() == 'stop':
            break
    
else:
    
    while True:
        accession = input("Enter an accession number (example: A5YKK6): ")
        accession = accession.replace("(", "").replace(")", "")
        if accession == 'stop':
            break
        accession = accession.upper()
        png = os.path.join(png_folder, f"{accession}.png")

        render_accession(accession)
        clear_shell()
        print_info(accession)
        display(Image(filename=png, width=set_width_png))
