# 08/08/2022
- Compare Eggnog6-PfamFams vs Eggnog6 OGs at Euk, Bact, Arq level
- Script for compare: OGs_vs_PfamFams.py


# 05/08/2022
- overlaps_none.emapper.hmm_hits.tgz has overlaping domains
- New script for cleaning domains : clean_clan_overlaping.py
- Pfam35 has 19632
- There are 18483 Domains in Eggnog6,
- 1149 Domains thar are not in Eggnog6 -> virus¿? Should check  

# 02/08/2022
- Carlos sent results Eggnog6 vs Pfam35: /scratch/carloslurm/eggnog6_pfam35/results/overlaps_none.emapper.hmm_hits.tgz
- emapper version 2.1.9: 

```
nohup hmm_mapper.py --mp_start_method forkserver 
 --hmm_maxhits 0 --cut_ga --usemem --hmm_maxseqlen 60000 
-clean_overlaps ${overlaps}  --cpu 60 --num_servers 30 
--qtype seq -i "${infile}"  --dbtype hmmdb -d "${pfamfn}"  
-o "${outdir}/${infn}"  > "logs/overlaps_${overlaps}.${infn}.out" 2>&1
```