RNAP_GENE_NAMES = [
    "fecI",
    "fliA",
    "rpoH",
    "rpoA",
    "rpoB",
    "rpoC",
    "rpoN",
    "rpoD",
    "rpoE",
    "rpoS",
    "rpoF",
    "rpoZ"]


def get_coding_genetic_target_len_d(component_name_str, genes_df):
    genetic_target_len_d = dict()

    # Multiple gene mutation
    if ',' in component_name_str:
        # Get genes
        gene_l = component_name_str.split(',')
        for genetic_target in gene_l:
            if genetic_target in NON_REGULONDB_GENE_L:
                component_len = get_non_regulonDB_gene_len(genetic_target)
            else:
                component_len = len(genes_df[genes_df["GENE_NAME"] == genetic_target]["GENE_SEQUENCE"].iloc[0])
            genetic_target_len_d[genetic_target] = component_len
        # Get intergenic region
        idx = 0
        while idx < len(gene_l):
            # building intergenic region annotation
            intergenic_region_str = gene_l[idx] + '/' + gene_l[idx+1]
            idx+=2
            d = get_intergenic_len_d(intergenic_region_str, genes_df)
            if d[intergenic_region_str] != 0:
                genetic_target_len_d.update(d)
    elif ';' in component_name_str:  # The case of pseudogenes
        pesudogene_l = component_name_str.split(';')
        for genetic_target in pesudogene_l:
            component_len = len(genes_df[genes_df["GENE_NAME"] == genetic_target]["GENE_SEQUENCE"].iloc[0])
            genetic_target_len_d[genetic_target] = component_len
    else:
        if component_name_str in NON_REGULONDB_GENE_L:
            component_len = get_non_regulonDB_gene_len(component_name_str)
        else:
            component_len = len(genes_df[genes_df["GENE_NAME"] == component_name_str]["GENE_SEQUENCE"].iloc[0])
        genetic_target_len_d[component_name_str] = component_len
        
    return genetic_target_len_d


def get_intergenic_len_d(component_name_str, genes_df):
    intergenic_len = 0
    gene_l = component_name_str.split('/')
    gene_1_pos_d = get_gene_pos_d(gene_l[0], genes_df)
    gene_2_pos_d = get_gene_pos_d(gene_l[1], genes_df)
    if gene_1_pos_d["right"] == gene_2_pos_d["left"] or gene_2_pos_d["right"] == gene_1_pos_d["left"]:  
        intergenic_len = 0  # when no intergenic region between 2 genes
    elif gene_l[0] == gene_l[1]:  # same gene for whatever reason; see with "gatC/gatC" in 42C exp.
        intergenic_len = len(genes_df[genes_df["GENE_NAME"] == gene_l[0]]["GENE_SEQUENCE"].iloc[0])
    elif gene_1_pos_d["right"] < gene_2_pos_d["left"]:
        intergenic_len = gene_2_pos_d["left"] - gene_1_pos_d["right"] - 1
    else:
        intergenic_len = gene_1_pos_d["left"] - gene_2_pos_d["right"] - 1
    return {component_name_str: intergenic_len}


NON_REGULONDB_GENE_L = ["ykfN", "insZ", "ydbA", "yehH", "yhiS", "yjgX", "ybeM", "ykgP"]


def get_non_regulonDB_gene_len(gene_name):
    gene_len = 0
    if gene_name == "ykfN":
        gene_len = 63  # https://biocyc.org/gene?orgid=ECOLI&id=G0-16715
    elif gene_name == "insZ":
        gene_len = 897  # https://ecocyc.org/gene?orgid=ECOLI&id=G6632
    elif gene_name == "ydbA":
        gene_len = 6063  # https://ecocyc.org/gene?orgid=ECOLI&id=G6632
    elif gene_name == "yehH":
        gene_len = 2860  # https://biocyc.org/gene?orgid=ECOLI&id=G8208
    elif gene_name == "yhiS":
        gene_len = 1224  # https://biocyc.org/gene?orgid=ECOLI&id=G8208
    elif gene_name == "yjgX":
        gene_len = 1199  # https://biocyc.org/gene?orgid=ECOLI&id=G7898
    elif gene_name == "ybeM":
        gene_len = 788  # https://biocyc.org/gene?orgid=ECOLI&id=G7898
    elif gene_name == "ykgP":
        gene_len = 90  # https://biocyc.org/gene?orgid=ECOLI&id=G7898
    return gene_len


def get_non_regulonDB_gene_pos(gene_name):
    gene_pos_t = (0,0)
    if gene_name == "ykfN":
        gene_pos_t = (263150, 263212)
    elif gene_name == "insZ":
        gene_pos_t = (1294426, 1295322)
    elif gene_name == "ydbA":
        gene_pos_t = (1465392, 1474013)
    elif gene_name == "yehH":
        gene_pos_t = (2197410, 2200269)
    elif gene_name == "yhiS":
        gene_pos_t = (3651291, 3653713)
    elif gene_name == "yjgX":
        gene_pos_t = (4499593, 4500791)
    elif gene_name == "ybeM":
        gene_pos_t = (65831, 658818)
    elif gene_name == "ykgP":
        gene_pos_t = (313716, 313805)
    return gene_pos_t
    
    
def get_gene_pos_d(gene_name, genes_df):
    gene_pos_d = dict()
    if gene_name in NON_REGULONDB_GENE_L:
        gene_pos_t = get_non_regulonDB_gene_pos(gene_name)
        gene_pos_d = {
            "left": gene_pos_t[0],
            "right": gene_pos_t[1]
        }
    else:
        gene_pos_d = {"left": genes_df[genes_df["GENE_NAME"]==gene_name]["GENE_POSLEFT"].iloc[0],
                      "right": genes_df[genes_df["GENE_NAME"]==gene_name]["GENE_POSRIGHT"].iloc[0]}
    return gene_pos_d


def get_gene_bnum(RegDB_obj_ID, regdb_gene_synonym_df):
    gene_bnum = ""
    g_synonym_df = regdb_gene_synonym_df[regdb_gene_synonym_df["OBJECT_ID"]==RegDB_obj_ID]
    g_bnum_df = g_synonym_df[g_synonym_df["OBJECT_SYNONYM_NAME"].str.contains('^b\d{4}')]
    if len(g_bnum_df) > 0:
        gene_bnum = g_bnum_df.iloc[0]["OBJECT_SYNONYM_NAME"]
    return gene_bnum
