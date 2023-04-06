#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Adriano De Marino, PhD
# Contact: adriano.demarino@gmail.com

import os 
from cyvcf2 import VCF
import pandas as pd
import argparse


def detect_build(vcf_file: str) -> str:
    """ Detect build of a VCF.
    Use the coordinates of common variants to identify the build / assembly of a genotype file
    that is being loaded.
    Notes
    -------
    References
    ----------
    1. Yates et. al. (doi:10.1093/bioinformatics/btu613),
       `<http://europepmc.org/search/?query=DOI:10.1093/bioinformatics/btu613>`_
    2. Zerbino et. al. (doi.org/10.1093/nar/gkx1098), https://doi.org/10.1093/nar/gkx1098
    3. Sherry ST, Ward MH, Kholodov M, Baker J, Phan L, Smigielski EM, Sirotkin K.
       dbSNP: the NCBI database of genetic variation. Nucleic Acids Res. 2001
       Jan 1;29(1):308-11.
    4. Database of Single Nucleotide Polymorphisms (dbSNP). Bethesda (MD): National Center
       for Biotechnology Information, National Library of Medicine. dbSNP accession: rs3094315,
       rs11928389, rs2500347, rs964481, rs2341354, rs3850290, and rs1329546
       (dbSNP Build ID: 151). Available from: http://www.ncbi.nlm.nih.gov/SNP/
    """
    cpus = os.cpu_count()
    rsids = [
        "rs3094315",
        "rs11928389",
        "rs2500347",
        "rs964481",
        "rs2341354",
        "rs3850290",
        "rs1329546",
    ]
    df = pd.DataFrame(
        {
            37: [
                'chr1_752566_G_A,T',
                'chr3_50927009_T_A,C,G',
                'chr1_144938320_T_A,C,G',
                'chrX_27656823_A_G,T',
                'chr1_918573_A_G',
                'chr14_23245301_T_A,C,G',
                'chrX_135474420_C_A,G,T',
            ],
            38: [
                'chr1_817186_G_A,T',
                'chr3_50889578_T_A,C,G',
                'chr1_148946169_T_A,C,G',
                'chrX_27638706_A_G,T',
                'chr1_983193_A_G',
                'chr14_22776092_T_A,C,G',
                'chrX_136392261_C_A,G,T',
            ],
        },
        index=rsids,
    )
    vcf = VCF(vcf_file, threads=cpus)
    build_keys = []
    ref_keys = df[37].tolist()+df[38].tolist()
    
    contig_name = vcf.raw_header.split('##contig=<ID=')[1].split('>\n')[0]
    
    ranges = [
        'chr1:752566-752566', 
        'chr3:50927009-50927009', 
        'chr1:144938320-144938320',
        'chrX:27656823-27656823',
        'chr1:918573-918573',
        'chr14:23245301-23245301',
        'chrX:135474420-135474420', 
        'chr1:817186-817186', 
        'chr3:50889578-50889578', 
        'chr1:148946169-148946169',
        'chrX:27638706-27638706',
        'chr1:983193-983193',
        'chr14:22776092-22776092',
        'chrX:136392261-136392261', 
    ]
    
    if 'chr' not in contig_name:
        ranges = [x.replace('chr','') for x in ranges]

    for ran in ranges:
        for variant in vcf(ran):
            target_chr = str(variant.CHROM)
            if "chr" not in target_chr:
                target_chr = "chr"+target_chr
            target_pos = str(variant.POS)
            target_ref = variant.REF
            target_alt = variant.ALT
            target_prefix = ("_").join([target_chr,target_pos,target_ref])
            for target_allele in target_alt:
                unique_id = target_prefix+"_"+target_allele
                for key in ref_keys:
                    ref_chr,ref_pos,ref_ref,ref_alt = (key).split("_")
                    ref_prefix = ("_").join([ref_chr,ref_pos,ref_ref])
                    ref_suffix = ref_alt.split(",")
                    for ref_alt in ref_suffix:
                        kk = ref_prefix+"_"+ref_alt
                        if kk == unique_id:
                            build_keys.append(key)
    vcf.close()

    if len(build_keys) == 0:
        print("No target variants found to detect the build")
        exit(1)

    elif len(df[df[37].isin(build_keys)]) >= 1 and len(df[df[38].isin(build_keys)]) == 0:
        return "Genome build: hg19/GRCh37"

    elif len(df[df[38].isin(build_keys)]) >= 1 and len(df[df[37].isin(build_keys)]) == 0:
        return "Genome build: hg38/GRCh38"

    else:
        print("Error in build detection")
        exit(1)


if __name__ == "__main__":

    argparser = argparse.ArgumentParser()
    argparser.add_argument('vcf_file', help='VCF file to detect the build')

    args = argparser.parse_args()
    vcf_file = args.vcf_file

    print(detect_build(vcf_file))
