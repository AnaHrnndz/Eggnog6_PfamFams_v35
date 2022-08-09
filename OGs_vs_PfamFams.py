#!/usr/bin/env python
from collections import defaultdict


dom2seqs = defaultdict(set)
total_seqs_pfam = set()
total_pfam_names = set()
with open('dom2seqs.tsv') as f:
    for line in f:
        if not line.startswith('#'):
            info = line.strip().split('\t')
            dom_name = info[0]
            seqs = set(info[4].split(','))
            dom2seqs[dom_name] = seqs
            total_seqs_pfam.update(seqs)
            total_pfam_names.add(dom_name)
print(len(dom2seqs))
print('SEQS in PFAM: ', len(total_seqs_pfam))

ogs2seqs = defaultdict(set)
total_seqs_og = set()
total_og_names = set()
with open('reduce_egg6_ogs_info_nogids.tsv') as f:
    for line in f:
        info = line.strip().split('\t')
        tax_level = info[0]
        og_name = info[2]
        seqs = set(info[6].split(','))
        ogs2seqs[og_name] = seqs
        total_seqs_og.update(seqs)
        total_og_names.add(og_name)
print(len(ogs2seqs))
print('SEQS in OG: ', len(total_seqs_og))



og_with_pfam = set()
pfam_with_og = set()

og_vs_pfam_out = open('og_vs_pfam.tsv', 'w')
pfam_vs_og_out = open('pfam_vs_og.tsv','w')
for og, og_seqs in ogs2seqs.items():
    for pfam, pfam_seqs in dom2seqs.items():
        common_seqs = list(og_seqs.intersection(pfam_seqs))
        if len(common_seqs) > 0:
            og_vs_pfam_out.write('\t'.join(map(str, [og, pfam, len(og_seqs), len(pfam_seqs), len(common_seqs), round(len(common_seqs)/len(og_seqs), 4), ','.join(common_seqs), '\n'])))
            pfam_vs_og_out.write('\t'.join(map(str, [pfam, og, len(pfam_seqs), len(og_seqs), len(common_seqs), round(len(common_seqs)/len(pfam_seqs), 4), ','.join(common_seqs), '\n'])))  

            og_with_pfam.add(og)
            pfam_with_og.add(pfam)
            

og_vs_pfam_out.close()
pfam_vs_og_out.close()

og_without_pfam = total_og_names.difference(og_with_pfam)
pfam_without_og = total_pfam_names.difference(pfam_with_og)

with open('og_without_pfam.tsv', 'w') as f:
    for line in og_without_pfam:
        f.write(line+'\n')

with open('pfam_without_og.tsv', 'w') as f:
    for line in pfam_without_og:
        f.write(line+'\n')