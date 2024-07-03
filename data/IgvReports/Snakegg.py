#!/usr/bin/env python3

# THIS SCRIPT SHOULD BE USED ONCE BEFORE THE SNAKEMAKE FILE
# FILL UP THE .env file first!

# IMPORTS ------------------------------------------------------------------

import os
import json
from dotenv import load_dotenv
from itertools import combinations
from collections import Counter
import pdb

load_dotenv() 

# FUNC ------------------------------------------------------------------

def get_all_pairs(fn):
    #====
    pairings=[]
    with open(fn) as f:
        for line in f:
            if line.startswith("sample"):
                header = line.rstrip().split(",")
                continue
            line = line.rstrip()
            fields = line.split(',')
            pairings.append(fields)

    # not scalable, lazy person! might actually be scalable though lol!
    full_pair = [j for j in pairings if j[1] != "" and j[2] != ""]
    ## need two sisters
    need_of_surr_gen1 = [i for i in pairings if i[1] == "" and i[2] == ""]
    ## need one sister
    need_of_surr_gen2 = [i for i in pairings if ((i[1] == "") ^ (i[2] == ""))]
    #===
    extra_pairs_list_gen2 = []
    for i in need_of_surr_gen2:
        all_same_gen = [j[0] for j in pairings if j[3] == i[3] and j[0] != i[0]]
        nmiss = sum([k == "" for k in i[1:3]])
        if nmiss == 1:
            for l in all_same_gen:
                #### ???
                idx_miss = [x for x in [1,2] if i[x] == ""][0]
                inst_add = [i[0], "", "", i[3]]
                inst_add[idx_miss] = l
                if idx_miss == 1:
                    inst_add[2] = i[2]
                else:
                    inst_add[1] = i[1]
                extra_pairs_list_gen2.append(inst_add)
        else:
            raise ValueError("jdbfsjhbfsj")
    #===
    all_comb = list(combinations([i[0] for i in need_of_surr_gen1], 2))
    extra_pairs_list_gen1 = []
    for i in need_of_surr_gen1:
        for j in all_comb:
            if i[0] in j:
                continue
            else:
                inst_add = [i[0], j[0], j[1], i[3]]
                extra_pairs_list_gen1.append(inst_add)
    all_pairings = full_pair + extra_pairs_list_gen1 + extra_pairs_list_gen2
    return all_pairings

def get_bam_files(path):
    abs_path = os.path.abspath(path)
    files = [i for i in os.listdir(abs_path) if not i.startswith('.')]
    files_bam = [i for i in files if i.endswith(".bam")]
    codes = [i.replace(".bam", "") for i in files_bam]
    codes2 = [i.split('.')[0] for i in codes]

    files_bam_dict = dict(zip(codes, files_bam))
    files_bai = [i for i in files if i.endswith(".bai")]
    codes_bai = [i.replace(".bai", "") for i in files_bai]
    files_bai_dict = dict(zip(codes_bai, files_bai))

    files_bam_ord = [files_bam_dict[i] for i in codes]
    files_bai_ord = [files_bai_dict[i] for i in codes]

    files_bam_abs = [os.path.join(abs_path, i) for i in files_bam_ord]
    files_bai_abs = [os.path.join(abs_path, i) for i in files_bai_ord]
    
    return codes, codes2, files_bam_abs, files_bai_abs

# files = C11.bed
def get_bed(path):
    abs_path = os.path.abspath(path)
    files = os.listdir(abs_path)
    samples = [os.path.splitext(i)[0] for i in files]
    files_abs = [os.path.join(abs_path, i) for i in files]
    return dict(zip(samples, files_abs))

### PARAMS ---------

qual = 20
path_bamfiles = os.getenv("BAMPATH")
path_beds_minimal =  os.getenv("BEDPATH_MINIMAL")
path_beds_moderate =  os.getenv("BEDPATH_MODERATE")
path_beds_stringent =  os.getenv("BEDPATH_STRINGENT")
fn_samples_ped = os.getenv("PEDPATH")
main_vcf = os.getenv("MAINVCF")

# HELPERS ---------------------------------------------

source_linked_bed_path = "TMP.nobackup/src/linkedbeds/"
source_linked_bam_path = "TMP.nobackup/src/linkedbams/"

if os.path.exists(source_linked_bed_path) or os.path.exists(source_linked_bam_path):
    raise ValueError("Need to remove linked directories in TMP.nobackup/")

os.makedirs(source_linked_bed_path)
os.makedirs(source_linked_bam_path)
os.makedirs("TMP.nobackup/src/reference")

### linking the files (BED)
bed_dict_min = get_bed(path_beds_minimal)
bed_dict_mod = get_bed(path_beds_moderate)
bed_dict_str = get_bed(path_beds_stringent)

samples, samples2, list_bams, bai = get_bam_files(path_bamfiles)
for i in range(len(samples)):
    os.symlink(list_bams[i], "{}/{}.bam".format(source_linked_bam_path,samples2[i]))
    os.symlink(bai[i], "{}/{}.bam.bai".format(source_linked_bam_path,samples2[i]))

numeric_samples = [i.split(".")[1] for i in samples]

## generate a equivalence files
with open("TMP.nobackup/src/reference/equivalence.csv", "w") as f:
    for i in range(len(samples)):
        f.write("{},{}\n".format(samples2[i], numeric_samples[i]))

### generating the final build-targets (aka results) =============
all_pairs = get_all_pairs(fn_samples_ped)

## 2P targets
all_kids = [i[0] for i in all_pairs]
all_kids_count = Counter(all_kids)

pairs_2p = [j for j in all_pairs if all_kids_count[j[0]] == 1]
pairs_3p = [j for j in all_pairs if all_kids_count[j[0]] == 2]
pairs_4p = [j for j in all_pairs if all_kids_count[j[0]] == 3]

target_codes_P2 = []
outpath_2P = "{}".format(source_linked_bed_path)
#os.makedirs(outpath_2P)
for i in pairs_2p:
    outfile_minimal = "{}_{}_{}_DNM_minimal.bed".format(i[0], i[1], i[2])
    outfile_moderate = "{}_{}_{}_DNM_moderate.bed".format(i[0], i[1], i[2])
    outfile_stringent = "{}_{}_{}_DNM_stringent.bed".format(i[0], i[1], i[2])
    target_codes_P2.append("{}_{}_{}".format(i[0], i[1], i[2]))
    
    olink = os.path.join(outpath_2P, outfile_minimal)
    os.symlink(
        bed_dict_min[i[0]],
        olink
    )

    olink = os.path.join(outpath_2P, outfile_moderate)
    os.symlink(
        bed_dict_mod[i[0]],
        olink
    )

    olink = os.path.join(outpath_2P, outfile_stringent)
    os.symlink(
        bed_dict_str[i[0]],
        olink
    )

## 3P targets
outpath_3P = "{}".format(source_linked_bed_path)
target_codes_P3 = []
kids = list(set([k[0] for k in pairs_3p]))
for j in kids:
    pts_j = []
    for m in pairs_3p:
        if m[0] == j:
            pts_j = pts_j + m[1:3]
    pts_list = list(set(pts_j))
    all_list = [j] + pts_list
    outfile_base = "_".join(all_list)
    target_codes_P3.append(outfile_base)
    outfile_minimal = "{}_DNM_minimal.bed".format(outfile_base)
    outfile_moderate = "{}_DNM_moderate.bed".format(outfile_base)
    outfile_stringent = "{}_DNM_stringent.bed".format(outfile_base)
    os.symlink(
        bed_dict_min[j],
        os.path.join(outpath_3P, outfile_minimal)
    )
    os.symlink(
        bed_dict_mod[j],
        os.path.join(outpath_3P, outfile_moderate)
    )
    os.symlink(
        bed_dict_str[j],
        os.path.join(outpath_3P, outfile_stringent)
    )

# 4P targets
outpath_4P = "{}".format(source_linked_bed_path)
target_codes_P4 = []
kids = list(set([k[0] for k in pairs_4p]))
for j in kids:
    pts_j = []
    for m in pairs_4p:
        if m[0] == j:
            pts_j = pts_j + m[1:3]
    pts_list = list(set(pts_j))
    all_list = [j] + pts_list
    outfile_base = "_".join(all_list)
    target_codes_P4.append(outfile_base)
    outfile_minimal = "{}_DNM_minimal.bed".format(outfile_base)
    outfile_moderate = "{}_DNM_moderate.bed".format(outfile_base)
    outfile_stringent = "{}_DNM_stringent.bed".format(outfile_base)
    os.symlink(
        bed_dict_min[j],
        os.path.join(outpath_4P, outfile_minimal)
    )
    os.symlink(
        bed_dict_mod[j],
        os.path.join(outpath_4P, outfile_moderate)
    )
    os.symlink(
        bed_dict_str[j],
        os.path.join(outpath_4P, outfile_stringent)
    )

# main vcf ------------------------------------------------------------------

os.symlink(main_vcf, "TMP.nobackup/src/reference/main.vcf.gz")
os.symlink("{}.tbi".format(main_vcf),
           "TMP.nobackup/src/reference/main.vcf.gz.tbi")

# results ------------------------------------------------------------------

results_P2_min = ["results.nobackup/{}_igvreport_minimal.html".format(i) for i in target_codes_P2]
results_P3_min = ["results.nobackup/{}_igvreport_minimal.html".format(i) for i in target_codes_P3]
results_P4_min = ["results.nobackup/{}_igvreport_minimal.html".format(i) for i in target_codes_P4]

results_P2_mod = ["results.nobackup/{}_igvreport_moderate.html".format(i) for i in target_codes_P2]
results_P3_mod = ["results.nobackup/{}_igvreport_moderate.html".format(i) for i in target_codes_P3]
results_P4_mod = ["results.nobackup/{}_igvreport_moderate.html".format(i) for i in target_codes_P4]

results_P2_str = ["results.nobackup/{}_igvreport_stringent.html".format(i) for i in target_codes_P2]
results_P3_str = ["results.nobackup/{}_igvreport_stringent.html".format(i) for i in target_codes_P3]
results_P4_str = ["results.nobackup/{}_igvreport_stringent.html".format(i) for i in target_codes_P4]

# output ------------------------------------------------------------------

## they are not really P4
config_list = {
    "qual": qual,
    "build-targets-p2": results_P2_min + results_P2_mod + results_P2_str,
    "build-targets-p3": results_P3_min + results_P3_mod + results_P3_str + results_P4_min + results_P4_mod + results_P4_str,
}

## write to json
with open("config.json", "w") as f:
    json.dump(config_list, f, indent=4)
