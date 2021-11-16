#!/usr/bin/env python

"""
Get download links for ensembl genome/gtf file by querying Ensembl website with BeautifulSoup.
The output will be a .nf file containing results are a groovy map.

Usage:
python get_ensembl_genome_gtf_download_links.py

@author: liuh, huk
"""

import urllib3
from bs4 import BeautifulSoup
import re
import json

def get_links(http, ftp_site, release, genome_download_links):
    ftp_site_fasta = ftp_site + release + "/fasta/"
    r = http.request('GET', ftp_site_fasta)
    soup = BeautifulSoup(r.data, 'html.parser')
    spcies_href = soup.find_all('a', href = re.compile(r'^[^.]'))

    species = [sp_h['href'].rstrip("/") for sp_h in spcies_href]


    for sp in species:
        genome_download_links.setdefault(sp, {})
        for k in ["gtf_md5sum", "genome_md5sum", "gtf", "genome"]:
            genome_download_links[sp].setdefault(k, None)

    species_link = [ftp_site_fasta + sp_h['href'] + "dna/" for sp_h in spcies_href]

    for i, sp_link in enumerate(species_link):
        req = http.request('GET', sp_link)
        soup = BeautifulSoup(req.data, 'html.parser')
        primary_seq_href = soup.find_all('a', href = re.compile(r'dna.primary_assembly.fa.gz'))
        toplevel_seq_href = soup.find_all('a', href = re.compile(r'dna.toplevel.fa.gz'))
        md5sum_href = soup.find_all('a', href = re.compile(r'CHECKSUMS'))

        if md5sum_href:
            genome_download_links[species[i]]["genome_md5sum"] = sp_link + md5sum_href[0]['href']
        if primary_seq_href:
            genome_download_links[species[i]]["genome"] = sp_link + primary_seq_href[0]['href']
        elif toplevel_seq_href:
            genome_download_links[species[i]]["genome"] = sp_link + toplevel_seq_href[0]['href']

    ## get gtf and gtf md5sum
    ftp_site_gtf = ftp_site + release + "/gtf/"
    r = http.request('GET', ftp_site_gtf)
    soup = BeautifulSoup(r.data, 'html.parser')
    spcies_href = soup.find_all('a', href = re.compile(r'^[^.]'))
    species_link = [ftp_site_gtf + sp_h['href']  for sp_h in spcies_href]

    species = [sp_h['href'].rstrip("/") for sp_h in spcies_href]
    for i, sp_link in enumerate(species_link):
        req = http.request('GET', sp_link)
        soup = BeautifulSoup(req.data, 'html.parser')
        gtf_href = soup.find_all('a', href = re.compile(r'\d+.gtf.gz'))
        md5sum_href = soup.find_all('a', href = re.compile(r'CHECKSUMS'))
        if md5sum_href:
            genome_download_links[species[i]]["gtf_md5sum"] = sp_link + md5sum_href[0]['href']
        if gtf_href:
            genome_download_links[species[i]]["gtf"] = sp_link + gtf_href[0]['href']
    return genome_download_links

if __name__ == "__main__":

    ## vertebrate_release, plant_release, fungi_release, metazoa_release can be passed as command-line
    ## parameters
    genome_download_links = {}
    http = urllib3.PoolManager()

    ## vertebrates
    ftp_site = "http://ftp.ensembl.org/pub/"
    vertebrate_release = "release-104"
    genome_download_links = get_links(http, ftp_site, vertebrate_release, genome_download_links)

    ## plants
    ftp_site = "http://ftp.ensemblgenomes.org/pub/plants/"
    plant_release = "release-51"
    genome_download_links = get_links(http, ftp_site, plant_release, genome_download_links)

    ## fungi
    ftp_site = "http://ftp.ensemblgenomes.org/pub/fungi/"
    fungi_release = "release-51"
    genome_download_links = get_links(http, ftp_site, fungi_release, genome_download_links)

    ## other metazoa
    ftp_site = "http://ftp.ensemblgenomes.org/pub/metazoa/"
    metazoa_release = "release-51"
    genome_download_links = get_links(http, ftp_site, metazoa_release, genome_download_links)

    ## Dump the dict as json:
    with open("../assets/genome_ensembl.json", "w") as outfile:
        json.dump(genome_download_links, outfile)

    ## Save genome names to file:
    with open("../modules/local/genome_ensembl.nf", "w") as f:
        f.write("// Map: ucsc name to links:\n")
        f.write("def get_genome_ensembl() {\n")
        f.write("\tdef list_genome_ensembl = [\n")
        for genome in genome_download_links:
            f.write('"' + genome + '",')
        f.write("\t]\n")
        f.write("return list_genome_ensembl\n")
        f.write("}")
