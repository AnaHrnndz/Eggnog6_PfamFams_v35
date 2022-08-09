#!/usr/bin/env python

from collections import defaultdict
import re


def checkIfDuplicates(listOfElems):
    ''' Check if given list contains any duplicates '''
    repeat_clans = list()    
    for elem in listOfElems:
        if listOfElems.count(elem) > 1:
            repeat_clans.append(elem)
    
    repeat_clans_updated = set([value for value in repeat_clans if value != 'None'])
    return repeat_clans_updated

def get_best_dom(hit2eval):
    
    for seq, doms in hit2eval.items():
        best_dom = defaultdict()

        for d in doms:
            start = int(seq2info[seq][d]['start'])
            end = int(seq2info[seq][d]['end'])
            eval_ = float(seq2info[seq][d]['evalue'])
            pos = set()
            for p in range( start, end):
                pos.add(p)

            if not best_dom:
                best_dom['name'] = d
                best_dom['pos'] = pos
                best_dom['eval'] = eval_

            elif len(pos.intersection(best_dom['pos'])) and eval_ < best_dom['eval']:
                best_dom['name'] = d
                best_dom['pos'] = pos
                best_dom['eval'] = eval_


    #print('BEST: ', best_dom['name'], best_dom['eval'])
    return best_dom


dom2clan = defaultdict()
with open('Pfam-A.clans.tsv') as f:
    for line in f:
        info = line.strip().split('\t')
        dom2clan[info[3]] = info[1]

seq2info = defaultdict(dict)
seq2clans =   defaultdict(list)
with open('overlaps_none.emapper.hmm_hits') as f:
    for line in f:
        if not line.startswith('#'):
            info = line.strip().split('\t')
            seq_name = info[0]
            clean_name1 = re.sub('\|', '_', seq_name)
            clean_name2 = re.sub('\:|\;|\*|,|\"|\[|\]|\(|\)|\/', '_', clean_name1)
            clean_name3 = re.sub(r"\\", '_', clean_name2)
            
            seq2info[clean_name3][info[1]] = dict()
            
            seq2info[clean_name3][info[1]]['evalue'] = info[2]
            seq2info[clean_name3][info[1]]['start'] = info[7]
            seq2info[clean_name3][info[1]]['end'] = info[8]
            seq2info[clean_name3][info[1]]['cov'] = info[9]
            if info[1] in dom2clan.keys():
                
                if dom2clan[info[1]] == '':
                    seq2info[clean_name3][info[1]]['clan'] = 'None'
                    seq2clans[clean_name3].append('None')
                else:
                    seq2info[clean_name3][info[1]]['clan'] = dom2clan[info[1]]
                    seq2clans[clean_name3].append(dom2clan[info[1]])
                              


# for seq, all_doms in seq2info.items():
    # for dom, info in all_doms.items():
        # print(seq, dom, info)

doms2seqs = defaultdict(set)
seq2info_no_overlaping = defaultdict(dict)
for seq, clans in seq2clans.items():
    repeat_clans = checkIfDuplicates(clans)
    
    if len(repeat_clans) > 0:
        hit2eval = defaultdict(list)

        for dom, info in seq2info[seq].items():
            if info['clan'] in repeat_clans:
                hit2eval[seq].append(dom)
            else:
                seq2info_no_overlaping[seq][dom] = info
                doms2seqs[dom].add(seq)
    
        best_dom = get_best_dom(hit2eval)
        name_best_dom = best_dom['name']
        info2add = seq2info[seq][name_best_dom]
        
        seq2info_no_overlaping[seq][name_best_dom] = info2add
        doms2seqs[name_best_dom].add(seq)


    else:
        seq2info_no_overlaping[seq] = seq2info[seq]
        for dom, info in seq2info[seq].items():
            doms2seqs[dom].add(seq)

# print('##########')
# for seq, all_doms in seq2info_no_overlaping.items():
    # for dom, info in all_doms.items():
        # print(seq, dom, info)

with open('seq2pfams.tsv', 'w') as f:
    f.write('\t'.join(['#Seq_name', 'Pfam', 'evalue', 'start', 'end', 'cov', 'clan', '\n']))
    for seq, doms in seq2info_no_overlaping.items():
        for dom, info in doms.items():
            f.write('\t'.join(map(str, [seq, dom, '\t'.join(info.values()), '\n'])))

with open('dom2seqs.tsv', 'w') as f:
    f.write('\t'.join(['#Pfam', 'Num_sp', 'Num_seqs', 'sp_list', 'seq_list', '\n']))
    for dom, seqs in doms2seqs.items():
        sp_set = set()
        for seq in seqs:
            sp_set.add(seq.split('.')[0])

        f.write('\t'.join(map(str, [dom, len(sp_set), len(seqs), ','.join(list(sp_set)), ','.join(list(seqs)), '\n'])))

#pfams2seqs: pfam, num_sp, num_seq, list_sp, list_seqs