#!/usr/bin/env python3

import os
from pyliftover import LiftOver
import pdb
lo = LiftOver('hg38', 'hg19')

import csv
from pyfaidx import Fasta
from dataclasses import dataclass
from tqdm import tqdm
from pysam import VariantFile
from Bio.Seq import Seq

from dotenv import load_dotenv
load_dotenv()
fasta_ref = os.getenv("FASTA")
fasta_polarized = os.getenv("PATH2ANCESTOR") + "/hs_ancestor.fa"

@dataclass
class Variant:
    chr: str
    pos: int
    ref: str
    alt: str
    label: str
    ref_tnp: str
    anc_tnp: str
    anc_ref: str

@dataclass
class Sample:
    name: str
    variants: list

    def get_num_variants(self):
        return len(self.variants)

ancestor_reference = Fasta(fasta_ref)
ancestor_polarized = Fasta(fasta_polarized)

csv_files = os.listdir("raw.data.nobackup/igvreports/")
vcf_files = os.listdir("raw.data.nobackup/unisamplevcf/")
vcf_files = [i for i in vcf_files if i.endswith("minimal.vcf.gz")]

names = [i.split("_")[0] for i in csv_files]
csv_dict = dict(zip(names, csv_files))

names2 = [i.split("_")[0] for i in vcf_files]
vcf_dict = dict(zip(names2, vcf_files))

samples = []

for sample_name in tqdm(csv_dict.keys()):
    csv_fn = "raw.data.nobackup/igvreports/{}".format(csv_dict[sample_name])
    vcf_fn = "raw.data.nobackup/unisamplevcf/{}".format(vcf_dict[sample_name])
    f = VariantFile(vcf_fn,'r')
    c42=Sample(name=sample_name, variants=[])
    # Replace 'path/to/file.csv' with the actual file path
    with open(csv_fn, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if row[0] == '':
                continue
            type_idx = [i for i in range(len(row)) if row[i] in ['red', 'yellow', 'green']]
            if len(type_idx) == 0:
                pdb.set_trace()
            type_idx = type_idx[0]
            type = row[type_idx]
            chr = row[0]
            pos = int(row[1])
            vcf_calls = []
            for rec in f.fetch(contig = chr, start = pos-1, end = pos):
                vcf_calls.append(rec)
            if len(vcf_calls) != 1:
                raise ValueError
            vcf_call = vcf_calls[0]
            ref = vcf_call.ref
            alt = vcf_call.alts[0]
            ref_tnp = ancestor_reference[chr][pos-2:pos+1]
            ### first need to liftover the variants
            new_coord = lo.convert_coordinate(chr, pos)
            if len(new_coord) != 1:
                anc_tnp = 'NA'
                anc_ref = 'NA'
            else:
                new_coord = new_coord[0]
                new_chr = new_coord[0]
                if new_chr not in ancestor_polarized.keys():
                    anc_tnp = 'NA'
                    anc_ref = 'NA'
                else:
                    new_pos = new_coord[1]
                    anc_tnp = ancestor_polarized[new_chr][new_pos-2:new_pos+1].seq
                    anc_ref = ancestor_polarized[new_chr][new_pos-1].seq
            c42.variants.append(Variant(chr=chr,
                                        pos=pos,
                                        label=type,
                                        ref = ref,
                                        alt = alt,
                                        ref_tnp = ref_tnp.seq,
                                        anc_tnp = anc_tnp,
                                        anc_ref = anc_ref))
    samples.append(c42)

# typeFam,mispol,ms
type_fam_dict = {
    'C11': 'trios',
    'C12': 'trios',
    'C31': 'trios',
    'C32': 'trios',
    'P1': 'surrogates',
    'P2': 'surrogates',
    'P3': 'surrogates',
    'P4': 'surrogates',
    'C41': 'surrogates',
    'C42': 'surrogates',
    'C21': 'surrogates',
    'C22': 'surrogates',
    'C23': 'surrogates'
}

def simplify_ms(ref, ctx, alt):
    if ref in ['C', 'T']:
        ms = "{ctx}>{alt}".format(ctx=ctx, alt=alt)
    else:
        rvcomp = str(Seq(ctx).reverse_complement())
        rvcomp_alt = str(Seq(alt).reverse_complement())
        ms = "{ctx}>{alt}".format(ctx=rvcomp, alt=rvcomp_alt)
    return ms

with open('myhped_ms6.csv', 'w', newline='', encoding='utf-8') as file_out:
    writer = csv.writer(file_out, lineterminator='\n')
    # Write the header row
    writer.writerow(['Sample', 'Chrom', 'Position', 'Ref', 'Alt', 'Type', 'Reference', 'Ancestor', 'AncestorRef', 'typeFam', 'mispol', 'ms'])
    for i in samples:
        for j in i.variants:
            ## determine if mispol
            if j.anc_ref == 'NA' or j.anc_ref == "N":
                mispol = "ancestorNA"
            elif j.anc_ref in [".", "-"]:
                mispol = "ancestorMissing"
            elif j.anc_ref.upper() == j.alt.upper():
                mispol = "mis-polarized"
            elif j.anc_ref.upper() == j.ref.upper():
                mispol = "correct"
            elif j.anc_ref.upper() in ["A", "C", "T", "G"]:
                mispol = "other"
            else:
                mispol = "dsfjhjskfh"
            ms = simplify_ms(j.ref, j.ref_tnp.upper(), j.alt)
            writer.writerow([i.name,
                             j.chr,
                             j.pos,
                             j.ref,
                             j.alt,
                             j.label,
                             j.ref_tnp,
                             j.anc_tnp,
                             j.anc_ref,
                             type_fam_dict[i.name],
                             mispol,
                             ms
                             ])
