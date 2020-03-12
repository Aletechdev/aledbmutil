from mut import get_inv_size, get_con_size, get_sub_size, get_del_size, get_ins_size, get_amp_size, is_coding_mut, get_MOB_type_str, get_codon_pos_chng, is_premature_stop_codon_SNP, is_start_codon_removal, get_SNP_aa_pos, is_disruptive, STRUCTURAL_LEVEL, get_original_nuc_mut_range


assert(is_coding_mut("intergenic (‑2/+1)")==False)
assert(is_coding_mut("A734V (GCG→GTG)")==True)
assert(is_coding_mut("noncoding (896/1542 nt)")==False)

assert(get_inv_size("1,443 bp inversion")==1443)

assert(get_con_size("3024 bp→REL606:592057‑588495")==3024)
assert(get_con_size("6 bp→REL606:1503176‑1503181")==6)

assert(get_sub_size("57511 bp→6 bp")==57505)
assert(get_sub_size('6827\xa0bp→81\xa0bp')==6746)
assert(get_sub_size('2\xa0bp→CT')==2)
assert(get_sub_size('2 bp→GA')==2)

assert(get_del_size("Δ82 bp") == 82)
assert(get_del_size("δ82 bp") == 82)
assert(get_del_size("Δ1,199 bp") == 1199)
assert(get_del_size("(T)60→50") == 10)

assert(get_ins_size("(TTC)1→2") == 3)
assert(get_ins_size("+GCTA") == 4)

assert(get_amp_size('36,044 bp x 3') == 72088)

assert(get_original_nuc_mut_range({"Mutation Type": "MOB", "Position": 99}) == (99,99))

assert(get_MOB_type_str("Δ1,199 bp")=="")
assert(get_MOB_type_str("IS1 (–) +8 bp")=="IS1")
assert(get_MOB_type_str("IS186 (–) +6 bp :: Δ1 bp")=="IS186")
assert(get_MOB_type_str("IS186 (–) +6 bp :: Δ1 bp")=="IS186")
assert(get_MOB_type_str("Δ1 :: IS186 (+) +6 bp :: Δ1")=="IS186")
test_str = "Δ6 bp :: REP161 (repetitive extragenic palindromic) element; contains 12 REP sequences (–) +2 bp :: +GGGGTGCCGCACTTCACAGCGGTGTAG"
assert(get_MOB_type_str(test_str)=="REP161")

assert(get_codon_pos_chng("ATC→AGC")==2)
assert(get_codon_pos_chng("AAC→AAT")==3)
assert(get_codon_pos_chng("AAC→AAC")==0)

assert(is_premature_stop_codon_SNP('P1100Q\xa0(CCG→TAG)\xa0'))

assert(get_SNP_aa_pos("P1100Q") == 1100)

assert(is_start_codon_removal('P1Q\xa0(CTG→CCG)\xa0'))
assert(is_start_codon_removal('P2Q\xa0(CTG→CCG)\xa0')==False)
assert(is_start_codon_removal('P1100Q\xa0(→TAG)\xa0')==False)
assert(is_start_codon_removal('M1M (ATG→ATA) †'))

# dictionaries can be substituted for Pandas series (DataFrame rows).
assert(is_disruptive({"coding": True, "genetic":True, "Mutation Type": "DEL", "target count": 1, "mutation size": 1}, STRUCTURAL_LEVEL)==True)
assert(is_disruptive({"coding": True, "genetic":True, "Mutation Type": "AMP", "target count": 1, "mutation size": 92}, STRUCTURAL_LEVEL)==True)

print("DONE")