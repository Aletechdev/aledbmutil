import pandas as pd
from gene import get_gene_bnum

df = pd.DataFrame.from_dict(
    {'OBJECT_ID': ['ECK120000001'],
     'OBJECT_SYNONYM_NAME': ['b4053'],
     'OS_INTERNAL_COMMENT': [None],
     'KEY_ID_ORG': ['ECK12']}, orient="columns")

assert(get_gene_bnum("ECK120000001", df) == "b4053")

print("DONE")