from mut import get_inv_size, get_con_size, get_sub_size, get_del_size, get_ins_size, get_amp_size, is_coding_mut, get_MOB_type_str, get_codon_pos_chng, is_premature_stop_codon_SNP, is_start_codon_removal, get_SNP_aa_pos, get_gene_count, get_clean_mut_gene_list, get_DEL_INS_MOB_aa_start_pos, get_DEL_AA_set, predict_mutation_effect_on_feature


i = "[araD],araA,araB"
e = ['araA', 'araB', 'araD']
assert(get_clean_mut_gene_list(i).sort() == e.sort())
assert(get_clean_mut_gene_list('[rph],[rph]') == ['rph', 'rph'])
i = {"Gene": "[araD],araA,araB", "Details": ''}
assert(get_gene_count(i) == 3)

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
assert(get_ins_size("(45 bp)1→2") == 45)

assert(get_amp_size('36,044 bp x 3') == 72088)

# assert(get_original_nuc_mut_range({"Mutation Type": "MOB", "Position": 99}) == (99,99))

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

assert(is_start_codon_removal('P2Q\xa0(CTG→CCG)\xa0')==False)
assert(is_start_codon_removal('P1100Q\xa0(→TAG)\xa0')==False)
assert(is_start_codon_removal('M1M (ATG→ATA) †'))

# dictionaries can be substituted for Pandas series (DataFrame rows).
# testing pseudogene INS, DEL, MOB logic
m = {"Mutation Type":"DEL", "coding":False, "mutation size":1}
f = {"genetic":True, "feature type":"pseudogene"}
assert(predict_mutation_effect_on_feature(m, f)=="other")

# start codon removal
m = {"Mutation Type":"SNP", "coding":True, "Details":"M1I (ATG→ATA) †"}
f = {"genetic":True, "feature type":"gene"}
assert(predict_mutation_effect_on_feature(m, f)=="truncation")

# premature stop codon
m = {"Mutation Type":"SNP", "coding":True, "Details":"P1100Q (CCG→TAG)"}
f = {"genetic":True, "feature type":"gene"}
assert(predict_mutation_effect_on_feature(m, f)=="truncation")

m = {"Mutation Type":"SNP", "coding":True, "Details":"A734V (GCG→GTG)"}
f = {"genetic":True, "feature type":"gene"}
assert(predict_mutation_effect_on_feature(m, f)=="nonsynonymous")

m = {"Mutation Type":"SNP", "coding":True, "Details":"A734A (GCG→GTG)"}
f = {"genetic":True, "feature type":"gene"}
assert(predict_mutation_effect_on_feature(m, f)=="synonymous")

for mt in ["AMP", "CNV", "SUB", "INV"]:
    assert(predict_mutation_effect_on_feature({"Mutation Type":"mt"}, {"feature type":"gene"})=="other")
    
m = {"Mutation Type":"INS", "mutation size":3}
f = {"genetic":True, "feature type":"gene"}
assert(predict_mutation_effect_on_feature(m, f)=="other")

m = {"Mutation Type":"INS", "mutation size":2}
f = {"genetic":True, "feature type":"gene"}
assert(predict_mutation_effect_on_feature(m, f)=="truncation")

m = {"Mutation Type":"INS", "mutation size":2}
f = {"genetic":False, "feature type":"promoter"}
assert(predict_mutation_effect_on_feature(m, f)=="other")

m = {"Mutation Type":"INS", "mutation size":10}
f = {"genetic":False, "feature type":"TFBS"}
assert(predict_mutation_effect_on_feature(m, f)=="truncation")

# testing output for "unknown" (intergenic region)
f = {"feature type":"unknown"}
assert(predict_mutation_effect_on_feature(None, f)=="other")

assert(get_DEL_INS_MOB_aa_start_pos("coding (1/20 nt)")==1)
assert(get_DEL_INS_MOB_aa_start_pos("coding (4/20 nt)")==2)
assert(get_DEL_INS_MOB_aa_start_pos("coding (58‑61/1413 nt)")==20)

# The below also tests get_DEL_AA_range
assert(get_DEL_AA_set("coding (1/20 nt)")=={1})  # test the removal of 1 BP from AA1
assert(get_DEL_AA_set("coding (2/20 nt)")=={1})  # test the removal of 1 BP from AA1
assert(get_DEL_AA_set("coding (3/20 nt)")=={1})  # test the removal of 1 BP from AA1
assert(get_DEL_AA_set("coding (4/20 nt)")=={2})
assert(get_DEL_AA_set("coding (1‑2/20 nt)")=={1})  # test removal of BPs within AA1
assert(get_DEL_AA_set("coding (1‑2/20 nt)")=={1})  # test removal of BPs within AA1
assert(get_DEL_AA_set("coding (1‑3/20 nt)")=={1})  # test removal of BPs within AA1
assert(get_DEL_AA_set("coding (4‑6/20 nt)")=={2})  # test of removal of all BPs for one AA
assert(get_DEL_AA_set("coding (4‑9/20 nt)")=={2,3})  # test for removal of all BPs for multiple AAs.
assert(get_DEL_AA_set("coding (6‑9/20 nt)")=={2,3})  # test for removal of some BPs for multiple AAs.
assert(get_DEL_AA_set("coding (6‑8/20 nt)")=={2,3})  # test for removal of some BPs for multiple AAs.
assert(get_DEL_AA_set("coding (4‑11/20 nt)")=={2,3,4})  # test for removal of some BPs for 3 AAs.
assert(get_DEL_AA_set("coding (2412‑2420/2547 nt")=={804, 805, 806, 807})  # 2412 is divisible by 3 == final nuc in an AA.
print("DONE")
