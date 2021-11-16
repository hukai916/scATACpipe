#!/usr/bin/env python

"""
Get download links for UCSC genome/gtf file by querying UCSC website with BeautifulSoup.
The output will be a .nf file containing results are a groovy map.

Usage:
python get_ucsc_genomes_with_gtf.py

@author: liuh, huk
"""

import urllib3
from bs4 import BeautifulSoup
import re
import copy
import json

lastest_download = "https://hgdownload.soe.ucsc.edu/downloads.html"
http = urllib3.PoolManager()
r = http.request('GET', lastest_download)
soup = BeautifulSoup(r.data, 'html.parser')

base="https://hgdownload.soe.ucsc.edu/"
bigzip_links = soup.find_all('a', href = re.compile("bigZips/$"))

genome_name = [link['href'].split("/")[2] for link in bigzip_links]
genome_fa_pattern = [re.compile(r'{}.fa.gz'.format(name)) for name in genome_name]

genome_download_links = {}
for gn in genome_name:
    genome_download_links.setdefault(gn, {})
    for k in ["md5sum", "gtf", "genome"]:
        genome_download_links[gn].setdefault(k, None)

for i, anchor in enumerate(bigzip_links):
    bgzp_link = base + anchor['href']

    #bgzp_link = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/"
    req = http.request('GET', bgzp_link)
    soup = BeautifulSoup(req.data, 'html.parser')
    md5sum_m = soup.find('a', href = re.compile(r'md5sum.txt'))
    gtf_m = soup.find('a', href = re.compile(r"^genes/$"))
    genome_m = soup.find('a', href = genome_fa_pattern[i])
    chrom_m = soup.find('a', href = re.compile(r"(chromFa.tar.gz|chromFa.zip)"))
    sm_chrom_m = soup.find('a', href = re.compile(r"softmask.fa.gz"))
    contig_scaffold_m = soup.find('a', href = re.compile(r"(contigFa.zip|scaffold.fa.gz)"))
    if md5sum_m:
        genome_download_links[genome_name[i]]["md5sum"] = bgzp_link + md5sum_m['href']
    if gtf_m:
        annotation_link = bgzp_link +gtf_m['href']
        gtf_res = http.request('GET', annotation_link)
        soup = BeautifulSoup(gtf_res.data, 'html.parser')
        ens_gtf = soup.find('a', href =re.compile("ensGene"))
        ref_gtf = soup.find('a', href =re.compile("refGene"))
        ncbiRef_gtf = soup.find('a', href =re.compile("ncbiRef"))
        ncbiGene_gtf = soup.find('a', href =re.compile("ncbiGenes"))

        ## if a genome has multiple gtf annotion files, the priority of gtf to be downloaded is
        ## ncbiRefSeq.gtf, ensGene.gtf, followed by refGene.gtf
        if ncbiRef_gtf:
            genome_download_links[genome_name[i]]["gtf"] = annotation_link + ncbiRef_gtf['href']
        elif ens_gtf:
            genome_download_links[genome_name[i]]["gtf"] = annotation_link+ ens_gtf['href']
        elif ref_gtf:
            genome_download_links[genome_name[i]]["gtf"] = annotation_link + ref_gtf['href']
        elif ncbiGene_gtf:
            genome_download_links[genome_name[i]]["gtf"] = annotation_link + ncbiGene_gtf['href']
    if genome_m:
        genome_download_links[genome_name[i]]["genome"] = bgzp_link + genome_m['href']
    elif chrom_m:
        genome_download_links[genome_name[i]]["genome"] = bgzp_link + chrom_m['href']
    elif sm_chrom_m:
        genome_download_links[genome_name[i]]["genome"] = bgzp_link + sm_chrom_m['href']
    elif contig_scaffold_m:
        genome_download_links[genome_name[i]]["genome"] = bgzp_link + contig_scaffold_m['href']

genome_download_links_copy = copy.deepcopy(genome_download_links)

# some species have no single genome.fa.gz file but a tar ball containing per chromosome fasta sequences
# some others have no genes/*.gtf annotation file
# remove genomes which have no a gtf annotation file and select the single genome.fa.gz file if possible,
#   otherwise choose the tar ball (Attention needed for furthur porcessing: uncompress first then concatenate)
for genome_name in genome_download_links_copy:
    if not genome_download_links_copy[genome_name]["gtf"] or not genome_download_links_copy[genome_name]["genome"]:
        del genome_download_links[genome_name]

genome_download_links   ## 146 genome assemblies

## Dump the dict as json:
with open("../assets/genome_ucsc.json", "w") as outfile:
    json.dump(genome_download_links, outfile)

# Output genome_ucsc.nf file that contains a list of genome names.
with open("../modules/local/genome_ucsc.nf", "w") as f:
    f.write("// Map: ucsc name to links:\n")
    f.write("def get_genome_ucsc() {\n")
    f.write("\tdef list_genome_ucsc = [\n")
    for genome in genome_download_links:
        f.write('"' + genome + '",')
    f.write("\t]\n")
    f.write("return list_genome_ucsc\n")
    f.write("}")
    #     hg38: [
    #       md5: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/md5sum.txt',
    #       gtf: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
    #       genome: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/hg38.fa.gz'
    #     ],
    #     hg19: [
    #       md5: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/md5sum.txt',
    #       gtf: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
    #       genome: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/hg38.fa.gz'
    #     ]
    #   ]
    #   return dict_genome_ucsc
    # }


    #     hg38: [
    #       md5: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/md5sum.txt',
    #       gtf: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
    #       genome: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/hg38.fa.gz'
    #     ],
    #     hg19: [
    #       md5: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/md5sum.txt',
    #       gtf: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz',
    #       genome: 'https://hgdownload.soe.ucsc.edu//goldenPath/hg38/bigZips/hg38.fa.gz'
    #     ]
    #   ]
    #   return dict_genome_ucsc
    # }




# for x in :
#     print(x, genome_download_links[x])
# final download links for each species
# with open("test4.txt", "w") as t_fh:
#     for g in sorted(genome_download_links):
#         t_fh.write(str(g)+ "\t")
#         for k in sorted(genome_download_links[g]):
#             t_fh.write(str(genome_download_links[g][k])+ "\t")
#         t_fh.write("\n")
