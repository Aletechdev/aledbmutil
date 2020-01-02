import re


# Could enhance code structure through strategy pattern: https://sourcemaking.com/design_patterns/strategy, https://sourcemaking.com/design_patterns/strategy/python/1


CODON_NUCLEOTIDE_COUNT = 3


def get_MOB_type_str(MOB_seq_change_str):
    MOB_type_str = ""
    MOB_seq_change_str = MOB_seq_change_str.replace(u'\xa0', u' ')
    MOB_str_ID = ""
    if "IS" in MOB_seq_change_str:
        MOB_str_ID = "IS"
    elif "REP" in MOB_seq_change_str:
        MOB_str_ID = "REP"
    
    if MOB_str_ID != "":
        IS_str_start_idx = MOB_seq_change_str.find(MOB_str_ID)
        IS_str_end_idx = MOB_seq_change_str.find(' ', IS_str_start_idx)
        MOB_type_str = MOB_seq_change_str[IS_str_start_idx:IS_str_end_idx]
        
    return MOB_type_str


def is_frameshift(nuc_shift_size):
    return nuc_shift_size % 3 != 0


def is_coding_mut(mut_details_str):
    is_coding = True
    noncoding_term_l = ["intergenic", "noncoding", "pseudogene"]
    if any(s in mut_details_str for s in noncoding_term_l):
        is_coding = False
    return is_coding


def is_genetic_mut(mut_details_str):
    is_genetic_mut = True
    nongenetic_term_l = ["intergenic"]
    if any(s in mut_details_str for s in nongenetic_term_l):
        is_genetic_mut = False
    return is_genetic_mut


def get_inv_size(seq_change_str):
    inv_size = 0
    seq_change_str = seq_change_str[:seq_change_str.find('bp')]
    inv_size = int(''.join([i for i in seq_change_str if i.isdigit()]))
    return inv_size


def get_con_size(seq_change_str):
    con_size = 0
    seq_change_str = seq_change_str[:seq_change_str.find('→')]
    con_size = int(''.join([i for i in seq_change_str if i.isdigit()]))
    return con_size


def get_MOB_size(seq_change_str):
    con_size = 0
    seq_change_str = seq_change_str[:seq_change_str.find('→')]
    con_size = int(''.join([i for i in seq_change_str if i.isdigit()]))
    return con_size


def get_sub_size(seq_change_str):
    sub_size = 0
    before_seq_change_str = seq_change_str[:seq_change_str.find('→')]
    before_seq_change_size = int(''.join([i for i in before_seq_change_str if i.isdigit()]))
    after_seq_change_str = seq_change_str[seq_change_str.find('→')+1:]
    if "bp" in after_seq_change_str:
        after_seq_change_size = int(''.join([i for i in after_seq_change_str if i.isdigit()]))
        sub_size = before_seq_change_size - after_seq_change_size  # expecting these types of substitutions to always shrink in size
    else:  # no change in size, but multiple base pairs replaced.
        sub_size = before_seq_change_size
    return sub_size


def _get_before_after_seq_change_size(seq_change_str):
    before_after_seq_change_size = 0
    if '→' in seq_change_str:
        s = seq_change_str[:seq_change_str.find('→')]
        before_size_str = ''.join([i for i in s if i.isdigit()])
        s = seq_change_str[seq_change_str.find('→')+1:]
        after_size_str = ''.join([i for i in s if i.isdigit()])
        before_after_seq_change_size = int(before_size_str) - int(after_size_str)
    return before_after_seq_change_size


def get_del_size(seq_change_str):
    del_size = 0
    if 'Δ' in seq_change_str or 'δ' in seq_change_str:
        del_size_str = ''.join([i for i in seq_change_str if i.isdigit()])
        del_size = int(del_size_str)
    if '→' in seq_change_str:
        del_size = _get_before_after_seq_change_size(seq_change_str)
    return del_size


def get_ins_size(seq_change_str):
    ins_size = 0
    if '→' in seq_change_str:
        s = seq_change_str[:seq_change_str.find('→')]
        before_seq_freq = int(''.join([i for i in s if i.isdigit()]))
        s = seq_change_str[seq_change_str.find('→')+1:]
        after_seq_freq = int(''.join([i for i in s if i.isdigit()]))
        if "bp" in seq_change_str:
            s = seq_change_str[seq_change_str.find('(')+1:seq_change_str.find(')')]
            s = s.replace(" bp", "")
            ins_size = int(s)
        else:
            seq_str = seq_change_str[seq_change_str.find('(')+1:seq_change_str.find(')')]
            ins_size = after_seq_freq * len(seq_str) - before_seq_freq * len(seq_str)
    if '+' in seq_change_str:
        ins_size = len(seq_change_str[seq_change_str.find('+')+1:])
    return ins_size


assert(get_ins_size("(TTC)1→2") == 3)
assert(get_ins_size("+GCTA") == 4)
assert(get_ins_size("(45 bp)1→2") == 45)


def get_amp_size(seq_change_str):
    seq_len_str = seq_change_str[:seq_change_str.find(' bp')]
    seq_len_str = seq_len_str.replace(',', '')
    seq_len = int(seq_len_str)
    multiplicity = int(seq_change_str[seq_change_str.find('x ')+len('x '):])
    seq_insertion_size = (multiplicity - 1) * seq_len
    return seq_insertion_size


def get_codon_pos_chng(codon_chng_str):
    ret_idx = 0
    codon_chng_list = codon_chng_str.split('→')
    for idx in range(CODON_NUCLEOTIDE_COUNT):  # assuming codon string is always going to be 3 long.
        if codon_chng_list[0][idx] != codon_chng_list[1][idx]:
            ret_idx = idx + 1
    return ret_idx


def is_in_TRN(details_str, gene_str, trn_gene_set):
    is_in_TRN = False
    if "intergenic" not in details_str:
        gene_list_str = gene_str
        gene_list_str = gene_list_str.replace('[','')
        gene_list_str = gene_list_str.replace(']','')
        gene_set = set(gene_list_str.split(', '))
        for gene in gene_set:
            if gene in trn_gene_set:
                is_in_TRN = True
    return is_in_TRN


# from https://en.wikipedia.org/wiki/Start_codon
def is_stop_codon(codon_str):
    codon_str = codon_str.lower()
    is_stop_codon = False
    stop_codon_l = ["taa", "tag", "tga"]
    if codon_str in stop_codon_l:
        is_stop_codon = True
    return is_stop_codon


# from https://en.wikipedia.org/wiki/Start_codon
def is_start_codon(codon_str):
    codon_str = codon_str.lower()
    is_start_codon = False
    start_codon_l = ["atg", "gtg", "ttg", "att", "ctg"]
    if codon_str in start_codon_l:
        is_start_codon = True
    return is_start_codon


def is_non_syn_SNP(amino_acid_change_str):
    return amino_acid_change_str[0] != amino_acid_change_str[-1]


def get_SNP_aa_pos(amino_acid_change_str):
    return int(re.sub("[^0-9]", "", amino_acid_change_str))


def get_codon_change_list(coding_SNP_details):
    codon_chng_str = coding_SNP_details[coding_SNP_details.find("(")+1:coding_SNP_details.find(")")]
    codon_change_list = codon_chng_str.split('→')
    return codon_change_list


def is_premature_stop_codon_SNP(coding_SNP_details):
    is_premature_stop_codon_SNP = False
    aa_chng_str = coding_SNP_details.split()[0]
    codon_chng_list = get_codon_change_list(coding_SNP_details)
    wt_codon = codon_chng_list[0]
    mut_codon = codon_chng_list[1]
    if is_non_syn_SNP(aa_chng_str) and is_stop_codon(mut_codon):
            is_premature_stop_codon_SNP = True
    return is_premature_stop_codon_SNP


def is_readthrough_codon_SNP(coding_SNP_details):
    is_readthrough_codon_SNP = False
    aa_chng_str = coding_SNP_details.split()[0]
    codon_change_list = get_codon_change_list(coding_SNP_details)
    wt_codon = codon_change_list[0]
    mut_codon = codon_change_list[1]
    if is_non_syn_SNP(aa_chng_str) \
        and is_stop_codon(wt_codon) and not is_stop_codon(mut_codon):
            is_readthrough_codon_SNP = True
    return is_readthrough_codon_SNP


def is_start_codon_removal(coding_SNP_details):
    is_start_codon_removal = False
    aa_chng_str = coding_SNP_details.split()[0]
    codon_change_list = get_codon_change_list(coding_SNP_details)
    wt_codon = codon_change_list[0]
    mut_codon = codon_change_list[1]
    if get_SNP_aa_pos(aa_chng_str) == 1 and is_start_codon(wt_codon) and not is_start_codon(mut_codon):
        is_start_codon_removal = True
    return is_start_codon_removal



regulatory_COG_description_set = {"Mobilome: prophages, transposons","Posttranslational modification, protein turnover, chaperones","Signal transduction mechanisms","Transcription","Translation, ribosomal structure and biogenesis"}

def has_regulatory_COG(cog_l_str):
    has_regulatory_COG = False
    cog_set = set(cog_l_str.split(';'))
    if (cog_set & regulatory_COG_description_set) != set():
        has_regulatory_COG = True
    return has_regulatory_COG


def is_regulatory(df_row):
    is_regulatory = False
    # Fast succeeding conditional structure
    if df_row["coding"] == False:
        is_regulatory = True
    elif df_row["TRN"] == True:
        is_regulatory = True
    else:
        cog_l_str = df_row["COG description"]
        is_regulatory = has_regulatory_COG(cog_l_str)
        
    return is_regulatory


def get_clean_mut_gene_list(gene_list_str):
    gene_list_str = gene_list_str.replace(" > ", " ")
    gene_list_str = gene_list_str.replace("]", "")
    gene_list_str = gene_list_str.replace("[", "")
    if "genes" in gene_list_str:
        start_idx = gene_list_str.rfind("genes")+len("genes")
        gene_list_str = gene_list_str[start_idx:]
    mut_gene_list = gene_list_str.split(", ")
    clean_mut_gene_list = [gene for gene in mut_gene_list]
    return clean_mut_gene_list


def get_gene_count(mut_df_row):
    target_count = 0
    gene_list_str = mut_df_row["mutation target annotation"]
    mut_details_str = mut_df_row["Details"]
    if "intergenic" in mut_details_str:
        target_count = 1
    else:
        gene_set = set(get_clean_mut_gene_list(gene_list_str))
        target_count = len(gene_set)
    return target_count



STRUCTURAL_LEVEL = 0
OPERATIONAL_LEVEL = 1


# TODO: should check to see if start codons are ever destroyed.
def is_disruptive_SNP(mut_df_row):
    is_disruptive_SNP = False
    if mut_df_row["Mutation Type"].lower() == "snp" and mut_df_row["coding"]:
        if is_premature_stop_codon_SNP(mut_df_row.Details) or is_readthrough_codon_SNP(mut_df_row.Details) or is_start_codon_removal(mut_df_row.Details):
            is_disruptive_SNP = True
    return is_disruptive_SNP


# TODO: Not currently checking if pseudogenes are being further disrupted. Pseudogenes can still be translates; need to check if they get further disrupted.
NON_CODING_COMP_COL_L = ['oriC', 'pseudogene', 'TFBS', 'promoter', 'TSS', 'RBS', 'cis-regulatory RNA', 'attenuator terminator', 'attenuator', 'terminator']
def is_disruptive(mut_df_row, metastruct_level):
    is_disruptive = False
    if metastruct_level == STRUCTURAL_LEVEL:
        if mut_df_row["genetic"]:
            if mut_df_row["Mutation Type"].lower() == "mob":
                is_disruptive = True
            # Assuming any type of mutation that affects multiple targets will remove a start or stop codon.
            elif mut_df_row["gene count"] > 1:
                is_disruptive = True
            elif mut_df_row["coding"]:
                if mut_df_row["Mutation Type"].lower() == "snp":
                    is_disruptive = is_disruptive_SNP(mut_df_row)
                elif any(x == mut_df_row["Mutation Type"].lower() for x in ["ins", "del", "amp"]):
                    is_disruptive = is_frameshift(mut_df_row["mutation size"])
        else:
            non_coding_comp_hit_l = mut_df_row[NON_CODING_COMP_COL_L].apply(lambda x: bool(x))
            if mut_df_row["Mutation Type"].lower() == "mob" and any(t for t in non_coding_comp_hit_l):
                is_disruptive = True
    else:
        print("!!! NEED TO IMPLEMENT !!!")
    
    return is_disruptive


# This is here and not in mut.py since it requires
# the parsing of multiple columns with a mutation row.
# Currently not returning the size of MOBs. Isn't something currently necessary.
def get_mut_size(mut_df_row):
    mut_size = 0  # Currently defaulting everything except for INS and DEL to 0 since don't need them.
    if mut_df_row["Mutation Type"] == "SNP":
        mut_size = 1
    elif mut_df_row["Mutation Type"] == "INS":
        mut_size = get_ins_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "DEL":
        mut_size = get_del_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "INV":
        mut_size = get_inv_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "CON":
        mut_size = get_inv_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "SUB":
        mut_size = get_sub_size(mut_df_row["Sequence Change"])
    elif mut_df_row["Mutation Type"] == "AMP":
        mut_size = get_amp_size(mut_df_row["Sequence Change"])
    return mut_size


# Returns the range of mutations to nucleotides in mutation region before mutation.
def get_original_nuc_mut_range(mut_df_row):
    mut_range = (0,0)
    if mut_df_row["Mutation Type"] == "SNP" \
    or mut_df_row["Mutation Type"] == "INS" \
    or mut_df_row["Mutation Type"] == "MOB" \
    or mut_df_row["Mutation Type"] == "AMP":
        mut_range = (mut_df_row["Position"], mut_df_row["Position"])
    elif mut_df_row["Mutation Type"] == "DEL" \
    or mut_df_row["Mutation Type"] == "INV" \
    or mut_df_row["Mutation Type"] == "CON" \
    or mut_df_row["Mutation Type"] == "SUB":
        mut_range = (mut_df_row["Position"], mut_df_row["Position"] - 1 + get_mut_size(mut_df_row))
    return mut_range
