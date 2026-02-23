#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 00:30:56 2024

@author: mechti
"""
import time
import requests
from concurrent.futures import ThreadPoolExecutor
import os, re, shutil, subprocess

# ----------------------------
# (B) TED rapid download: aria2 + validate + retry (only used for TED)
# ----------------------------
def download_structures_ted_rapid(ted_ids, directory, cores_num,
                                 base_url="https://ted.cathdb.info/api/v1/files"):
    os.makedirs(directory, exist_ok=True)
    if not ted_ids:
        return
    ids = ted_ids
    # Clean + dedupe preserving order
    seen = set()
    ids = []
    for x in ted_ids:
        x = str(x).rstrip("\r\n")
        x = re.sub(r"#.*$", "", x).strip()
        if not x:
            continue
        if x not in seen:
            ids.append(x)
            seen.add(x)
    if not ids:
        return

    if shutil.which("aria2c") is None:
        # keep AF/ESM unaffected; only TED requires aria2 for rapid mode
        raise RuntimeError("aria2c not found. Install aria2 (e.g., apt-get install aria2) to enable TED rapid download.")

    parent_dir = os.path.dirname(os.path.abspath(directory))
    work_dir = os.path.join(parent_dir, "ted_dl_work")  # sibling, not inside query folder
    os.makedirs(work_dir, exist_ok=True)

    urls_file = os.path.join(work_dir, "ted_urls.txt")
    retry_urls_file = os.path.join(work_dir, "retry_urls.txt")
    aria_log = os.path.join(work_dir, "aria2.log")
    missing_ids_file = os.path.join(work_dir, "missing_ids.txt")
    bad_files_file = os.path.join(work_dir, "bad_files.txt")

    def looks_bad(fp):
        if (not os.path.exists(fp)) or os.path.getsize(fp) == 0:
            return True
        with open(fp, "rb") as f:
            head = f.read(512).decode("utf-8", errors="ignore").lower()
        return any(m in head for m in ("<html", "<!doctype", "error", "not found", "too many requests", "forbidden"))

    # Build URL list
    with open(urls_file, "w", encoding="utf-8") as f:
        for tid in ids:
            f.write(f"{base_url}/{tid}.pdb\n")

    # Pass 1 (fast)
    jobs = max(8, int(cores_num) * 8)
    subprocess.run([
        "aria2c",
        f"--input-file={urls_file}",
        f"--dir={directory}",
        "-j", str(jobs),
        "--continue=true",
        "--max-connection-per-server=8",
        "--split=8",
        "--min-split-size=1M",
        "--file-allocation=none",
        "--retry-wait=3",
        "--max-tries=20",
        "--timeout=60",
        "--connect-timeout=20",
        "--auto-file-renaming=false",
        "--allow-overwrite=true",
        "--summary-interval=10",
        "--console-log-level=notice",
        f"--log={aria_log}",
    ], check=False)

    # Verify missing + bad
    missing = []
    bad = []
    for tid in ids:
        fp = os.path.join(directory, f"{tid}.pdb")
        if (not os.path.exists(fp)) or os.path.getsize(fp) == 0:
            missing.append(tid)
        elif looks_bad(fp):
            bad.append(tid)

    with open(missing_ids_file, "w", encoding="utf-8") as f:
        f.write("\n".join(missing) + ("\n" if missing else ""))
    with open(bad_files_file, "w", encoding="utf-8") as f:
        f.write("\n".join([f"{x}.pdb" for x in bad]) + ("\n" if bad else ""))

    # Retry missing+bad only (safer)
    if missing or bad:
        retry_urls = []
        retry_urls += [f"{base_url}/{tid}.pdb" for tid in missing]
        retry_urls += [f"{base_url}/{tid}.pdb" for tid in bad]
        retry_urls = list(dict.fromkeys(retry_urls))

        with open(retry_urls_file, "w", encoding="utf-8") as f:
            f.write("\n".join(retry_urls) + "\n")

        subprocess.run([
            "aria2c",
            f"--input-file={retry_urls_file}",
            f"--dir={directory}",
            "-j", str(max(4, jobs // 2)),
            "--continue=true",
            "--max-connection-per-server=4",
            "--split=4",
            "--min-split-size=1M",
            "--file-allocation=none",
            "--retry-wait=5",
            "--max-tries=30",
            "--timeout=60",
            "--connect-timeout=20",
            "--auto-file-renaming=false",
            "--allow-overwrite=true",
            "--summary-interval=10",
            "--console-log-level=notice",
            f"--log={aria_log}",
        ], check=False)

    present = sum(
        1 for tid in ids
        if os.path.exists(os.path.join(directory, f"{tid}.pdb"))
        and not looks_bad(os.path.join(directory, f"{tid}.pdb"))
    )
    print(f"[TED aria2] present {present}/{len(ids)} | log: {aria_log}")

# Generalized download function with file format handling
def download_files(url_template, ids, directory, file_extension):
    start_time=time.time()

    def download_file(file_id):
        filename = os.path.join(directory, f'{file_id}.{file_extension}')
        if not os.path.isfile(filename):
            url = url_template.format(file_id, file_extension)
            response = requests.get(url)
            if response.status_code == 200:
                with open(filename, 'wb') as file:
                    file.write(response.content)
                print(f'Downloaded {filename}')
            else:
                print(f'Failed to download {file_id}.{file_extension}')

    with ThreadPoolExecutor(max_workers=3) as executor:
        executor.map(download_file, ids)

    print("Overall Running Time:", time.time() - start_time)


def download_structures_af_rapid(uniprot_accessions, directory, cores_num,
                                base_url="https://alphafold.ebi.ac.uk/files",
                                model_version="v6",
                                try_gz_fallback=True):
    """
    Rapid AlphaFold download via aria2c.
    Downloads mmCIF files for UniProt accessions as:
      AF-<ACC>-F1-model_<model_version>.cif
    Optionally retries with .cif.gz if .cif is missing.
    """
    os.makedirs(directory, exist_ok=True)
    if not uniprot_accessions:
        return

    # Clean + dedupe preserving order
    seen = set()
    accs = []
    for x in uniprot_accessions:
        x = str(x).rstrip("\r\n")
        x = re.sub(r"#.*$", "", x).strip()
        if not x:
            continue
        if x not in seen:
            accs.append(x)
            seen.add(x)
    if not accs:
        return

    if shutil.which("aria2c") is None:
        raise RuntimeError("aria2c not found. Install aria2 (apt-get install aria2) to enable AF rapid download.")

    parent_dir = os.path.dirname(os.path.abspath(directory))
    work_dir = os.path.join(parent_dir, "af_dl_work")
    os.makedirs(work_dir, exist_ok=True)

    urls_file = os.path.join(work_dir, "af_urls.txt")
    retry_urls_file = os.path.join(work_dir, "retry_urls.txt")
    aria_log = os.path.join(work_dir, "aria2.log")
    missing_ids_file = os.path.join(work_dir, "missing_ids.txt")
    bad_files_file = os.path.join(work_dir, "bad_files.txt")

    def looks_bad(fp):
        # Detect HTML error pages / empty files etc.
        if (not os.path.exists(fp)) or os.path.getsize(fp) == 0:
            return True
        with open(fp, "rb") as f:
            head = f.read(512).decode("utf-8", errors="ignore").lower()
        return any(m in head for m in ("<html", "<!doctype", "error", "not found", "too many requests", "forbidden"))

    def af_fname(acc, ext="cif"):
        return f"AF-{acc}-F1-model_{model_version}.{ext}"

    # Build URL list (primary: .cif)
    with open(urls_file, "w", encoding="utf-8") as f:
        for acc in accs:
            f.write(f"{base_url}/{af_fname(acc, 'cif')}\n")

    # Pass 1
    # Be a bit gentler than TED (single host, rate limits can happen)
    jobs = max(8, int(cores_num) * 4)
    subprocess.run([
        "aria2c",
        f"--input-file={urls_file}",
        f"--dir={directory}",
        "-j", str(jobs),
        "--continue=true",
        "--max-connection-per-server=4",
        "--split=4",
        "--min-split-size=1M",
        "--file-allocation=none",
        "--retry-wait=3",
        "--max-tries=20",
        "--timeout=60",
        "--connect-timeout=20",
        "--auto-file-renaming=false",
        "--allow-overwrite=true",
        "--summary-interval=10",
        "--console-log-level=notice",
        f"--log={aria_log}",
    ], check=False)

    # Verify missing + bad
    missing = []
    bad = []
    for acc in accs:
        fp = os.path.join(directory, af_fname(acc, "cif"))
        if (not os.path.exists(fp)) or os.path.getsize(fp) == 0:
            missing.append(acc)
        elif looks_bad(fp):
            bad.append(acc)

    with open(missing_ids_file, "w", encoding="utf-8") as f:
        f.write("\n".join(missing) + ("\n" if missing else ""))
    with open(bad_files_file, "w", encoding="utf-8") as f:
        f.write("\n".join([af_fname(x, "cif") for x in bad]) + ("\n" if bad else ""))

    # Retry missing+bad (.cif)
    if missing or bad:
        retry_urls = []
        retry_urls += [f"{base_url}/{af_fname(acc, 'cif')}" for acc in missing]
        retry_urls += [f"{base_url}/{af_fname(acc, 'cif')}" for acc in bad]
        retry_urls = list(dict.fromkeys(retry_urls))

        with open(retry_urls_file, "w", encoding="utf-8") as f:
            f.write("\n".join(retry_urls) + "\n")

        subprocess.run([
            "aria2c",
            f"--input-file={retry_urls_file}",
            f"--dir={directory}",
            "-j", str(max(4, jobs // 2)),
            "--continue=true",
            "--max-connection-per-server=2",
            "--split=2",
            "--min-split-size=1M",
            "--file-allocation=none",
            "--retry-wait=5",
            "--max-tries=30",
            "--timeout=60",
            "--connect-timeout=20",
            "--auto-file-renaming=false",
            "--allow-overwrite=true",
            "--summary-interval=10",
            "--console-log-level=notice",
            f"--log={aria_log}",
        ], check=False)

    # Optional: fallback to .cif.gz for remaining missing
    if try_gz_fallback:
        still_missing = []
        for acc in accs:
            fp = os.path.join(directory, af_fname(acc, "cif"))
            if (not os.path.exists(fp)) or looks_bad(fp):
                still_missing.append(acc)

        if still_missing:
            gz_urls_file = os.path.join(work_dir, "af_urls_gz.txt")
            with open(gz_urls_file, "w", encoding="utf-8") as f:
                for acc in still_missing:
                    f.write(f"{base_url}/{af_fname(acc, 'cif.gz')}\n")

            subprocess.run([
                "aria2c",
                f"--input-file={gz_urls_file}",
                f"--dir={directory}",
                "-j", str(max(4, jobs // 2)),
                "--continue=true",
                "--max-connection-per-server=2",
                "--split=2",
                "--min-split-size=1M",
                "--file-allocation=none",
                "--retry-wait=5",
                "--max-tries=30",
                "--timeout=60",
                "--connect-timeout=20",
                "--auto-file-renaming=false",
                "--allow-overwrite=true",
                "--summary-interval=10",
                "--console-log-level=notice",
                f"--log={aria_log}",
            ], check=False)

    present = 0
    for acc in accs:
        fp_cif = os.path.join(directory, af_fname(acc, "cif"))
        fp_cifgz = os.path.join(directory, af_fname(acc, "cif.gz"))
        ok = (os.path.exists(fp_cif) and not looks_bad(fp_cif)) or (os.path.exists(fp_cifgz) and not looks_bad(fp_cifgz))
        present += 1 if ok else 0

    print(f"[AF aria2] present {present}/{len(accs)} | log: {aria_log}")


# Specific download functions for AF, PDB (in CIF), and esm (in PDB)
def download_structures_af(uniprot_accessions, directory, cores_num):
    # af_ids = ["AF-" + accession + "-F1-model_v6" for accession in uniprot_accessions]
    # url_template = 'https://alphafold.ebi.ac.uk/files/{}.{}'
    #
    # download_files(url_template, af_ids, directory, file_extension="cif")
    download_structures_af_rapid(uniprot_accessions, directory, cores_num)

def download_structures_pdb(pdb_ids, directory):
    url_template = 'https://files.rcsb.org/download/{}.{}'
    download_files(url_template, pdb_ids, directory, file_extension="cif")


def download_structures_esm(esm_ids, directory):
    url_template = 'https://api.esmatlas.com/fetchPredictedStructure/{}.{}'  # Replace with actual URL for esm PDBs
    download_files(url_template, esm_ids, directory, file_extension="pdb")

def download_structures_ted(ted_ids, directory,cores_num):
    # url_template='https://ted.cathdb.info/api/v1/files/{}.{}'
    # download_files(url_template, ted_ids, directory, file_extension="pdb")
    download_structures_ted_rapid(ted_ids, directory, cores_num)     # ✅ TED only uses rapid aria2 method; AF/ESM/PDB unaffected

def main(uniprot_accessions, pdb_ids, esm_ids, ted_ids, path_query_structures, cores_num):
    # AlphaFold and PDB downloads use CIF format
    download_structures_af(uniprot_accessions, path_query_structures,cores_num)
    download_structures_pdb(pdb_ids, path_query_structures)

    # esm and ted downloads use PDB format
    download_structures_esm(esm_ids, path_query_structures)
    download_structures_ted(ted_ids, path_query_structures, cores_num)



if __name__ == "__main__":
    cores_num=2
    uniprot_accessions = ['A0A2K5XT84', 'G3QSU8', 'A5Z1T7']
    pdb_ids = ['1CRN', '4HHB', '5XNL']
    esm_ids = ['MGYP002537940442', 'MGYP001215146166', 'MGYP001823580159']
    ted_ids=['AF-A0A002-F1-model_v4_TED01','AF-A0A003-F1-model_v4_TED01','AF-A0A009EAK7-F1-model_v4_TED01']
    main(uniprot_accessions, pdb_ids, esm_ids, ted_ids, path_query_structures=".", cores_num=cores_num)
