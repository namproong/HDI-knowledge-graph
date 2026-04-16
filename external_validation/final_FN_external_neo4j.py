#%%
import pandas as pd
import pickle
from neo4j import GraphDatabase
# %%
uri = "neo4j://127.0.0.1:7687"
username = "neo4j"
password = "HDIKGN2026"
driver = GraphDatabase.driver(uri, auth=(username, password))
print("Connected successfully!")
#%%
FN_compound = pd.read_csv("exFN_compound_level_evidence.csv")
print(len(FN_compound)) #182
print(FN_compound.head()) 
#%%
print(FN_compound.columns) #Index(['species_tax_id', 'species_name', 'compound_key_1', 'compound_name_1', 'compound_key_2', 'compound_name_2', 'probability', 'prediction'],
#%%
FN_compound = FN_compound.reset_index(drop=True)
FN_compound["pair_id"] = FN_compound.index
# %%
FN_pairs = list(zip(FN_compound.compound_key_1, FN_compound.compound_key_2))
print("FN_pairs:", len(FN_pairs)) #FN_pairs: 182
#%%
FN_pairs = FN_compound[
    ["pair_id","compound_key_1","compound_key_2","probability"]
].drop_duplicates()
#%%
num_unique_pairs = FN_compound[
    ["species_name", "compound_key_2"]
].drop_duplicates().shape[0]

print(num_unique_pairs) #6
#%% load neo4j node name 
hdi_target = pd.read_csv("hdi_target.csv")
print(hdi_target.columns) #Index(['target_key', 'source', 'source_target_id', 'tid', 'pref_name', 'target_type', 'organism', 'tax_id', 'is_protein'],      dtype='object')
compound_lookup = pd.read_csv("compound_lookup_master.csv")
compound_lookup = compound_lookup[[
    "compound_key",
    "display_name",
]]
target_lookup = hdi_target[["target_key","pref_name"]].rename(
    columns={
        "target_key": "target_id",
        "pref_name": "target_name"
    }
)
target_lookup = target_lookup.drop_duplicates(subset=["target_id"])
#%%
def run_query(tx, query, c1, c2):
    result = tx.run(query, c1=c1, c2=c2)
    return [dict(record) for record in result] 
# %%
#Query สำหรับ Target
TARGET_QUERY = """
MATCH (c1:Compound {id:$c1})-[r1]->(t:Target)<-[r2]-(c2:Compound {id:$c2})
RETURN DISTINCT
"Target" AS node_type,
t.id AS node_id,
type(r1) AS rel_c1,
type(r2) AS rel_c2
"""
#ATC hierarchy
ATC_QUERY = """
MATCH (c1:Compound {id:$c1})-[:HAS_ATC5]->(a5_1:ATC5)
MATCH (c2:Compound {id:$c2})-[:HAS_ATC5]->(a5_2:ATC5)

MATCH (a5_1)-[:SUBCLASS_OF_ATC4]->(a4:ATC4)
MATCH (a5_2)-[:SUBCLASS_OF_ATC4]->(a4:ATC4)

RETURN DISTINCT
"ATC" AS node_type,
a4.id AS node_id,
"SHARES_ATC" AS rel_c1,
"SHARES_ATC" AS rel_c2
"""
#Metabolite (compound → compound)
METABOLITE_QUERY = """
MATCH (c1:Compound {id:$c1})-[:METABOLIZED_TO]->(c2:Compound {id:$c2})
RETURN
"Metabolite_direct" AS node_type,
c2.id AS node_id,
"METABOLIZED_TO" AS rel_c1,
"METABOLIZED_FROM" AS rel_c2

UNION

MATCH (c2:Compound {id:$c2})-[:METABOLIZED_TO]->(c1:Compound {id:$c1})
RETURN
"Metabolite_direct" AS node_type,
c1.id AS node_id,
"METABOLIZED_FROM" AS rel_c1,
"METABOLIZED_TO" AS rel_c2
"""
# %%
all_queries = [
    TARGET_QUERY,
    ATC_QUERY,
    METABOLITE_QUERY
]

FN_evidence_rows = []

with driver.session(database="hdikgv3") as session:

    for _, row in FN_pairs.iterrows():
        pair_id = row["pair_id"]
        c1 = row["compound_key_1"]
        c2 = row["compound_key_2"]
        prob = row["probability"]

        for query in all_queries:

            results = session.execute_read(run_query, query, c1, c2)

            for row_q in results:
                FN_evidence_rows.append([
                    pair_id,
                    c1,
                    c2,
                    prob,
                    row_q["node_type"],
                    row_q["node_id"],
                    row_q["rel_c1"],
                    row_q["rel_c2"]
                ])
#%%
def classify_relation(rel, node_type):

    if node_type == "ATC":
        return "Drug class similarity"

    if node_type == "Metabolite_direct":
        return "Metabolic relationship"

    PD = [
        "NEGATIVELY_MODULATES",
        "POSITIVELY_MODULATES",
        "MODULATES",
        "BINDS",
        "AFFECTS_ENZYMATIC_FUNCTION",
        "STRUCTURALLY_MODULATES",
        "INTERACTS_WITH"
    ]

    ENZYME = [
        "INHIBITS",
        "INDUCES",
        "ACTIVATES",
        "DOWNREGULATES"
    ]

    PK = [
        "METABOLIZED_TO",
        "METABOLIZED_BY",
        "IS_TRANSPORT_SUBSTRATE_OF"
    ]

    if rel in PD:
        return "Pharmacodynamic interaction"
    elif rel in ENZYME:
        return "Enzyme regulation"
    elif rel in PK:
        return "Drug metabolism (PK)"
    else:
        return "Unclassified"
#%%
final_rows = []

for row in FN_evidence_rows:

    node_type = row[4]   
    rel1 = row[6]        
    rel2 = row[7] 

    final_rows.append(row + [
        classify_relation(rel1, node_type),
        classify_relation(rel2, node_type)
    ])
#%%
columns = [
    "pair_id",
    "herb_compound_key",
    "drug_compound_key",
    "probability",
    "node_type",
    "node_id",
    "rel_c1",
    "rel_c2",
    "mechanism_c1",
    "mechanism_c2" 
]

df_fp = pd.DataFrame(final_rows, columns=columns)
#%%
df_fp["node_id"] = df_fp["node_id"].astype(str)
target_lookup["target_id"] = target_lookup["target_id"].astype(str)
#%%
#%%
# merge ชื่อ herb 
df_fp = df_fp.merge(
    FN_compound[[
        "pair_id",
        "species_name"
    ]],
    on="pair_id",
    how="left"
)#%%
df_fp.head(5)
#%%
df_fp = df_fp.sort_values(
    by=["species_name","drug_compound_key","probability"],
    ascending=[True, True, False]
)
#%%
df_fp.head()
#%%
cols = [
    "pair_id",
    "species_name",
    "herb_compound_key",
    "drug_compound_key",
    "probability",
    "node_type",
    "node_id",
    "rel_c1",
    "rel_c2",
    "mechanism_c1",
    "mechanism_c2"
]

df_fp = df_fp[cols]
#%%
df_fp.to_csv("final_exFN_evidence_table.csv", index=False)
# %%
PD_REL = {
    "NEGATIVELY_MODULATES",
    "POSITIVELY_MODULATES",
    "MODULATES",
    "BINDS",
    "AFFECTS_ENZYMATIC_FUNCTION",
    "STRUCTURALLY_MODULATES",
    "INTERACTS_WITH"
}

ENZYME_REG = {
    "INHIBITS",
    "INDUCES",
    "ACTIVATES",
    "DOWNREGULATES"
}

PK_CORE = {
    "METABOLIZED_BY",
    "IS_TRANSPORT_SUBSTRATE_OF"
}
# %%
def derive_final_mechanism(row):

    r1 = row["rel_c1"]
    r2 = row["rel_c2"]
    node = row["node_type"]
    pair = {r1, r2}  

    if node == "Target":

        # --- PK: competition ---
        if r1 == "METABOLIZED_BY" and r2 == "METABOLIZED_BY":
            return "PK_competition"

        # --- PK: enzyme regulation ---
        if (pair & ENZYME_REG) and "METABOLIZED_BY" in pair:
            return "PK_enzyme_regulation"

        # --- PK: transporter ---
        if r1 == "IS_TRANSPORT_SUBSTRATE_OF" and r2 == "IS_TRANSPORT_SUBSTRATE_OF":
            return "PK_transporter"

        # --- PK: transporter inhibition ---
        if (pair & ENZYME_REG) and "IS_TRANSPORT_SUBSTRATE_OF" in pair:
            return "PK_transporter_inhibition"

        # --- PD: same target ---
        if r1 in PD_REL and r2 in PD_REL:
            return "PD_target"

        # --- enzyme modulation  
        if r1 in ENZYME_REG and r2 in ENZYME_REG:
            return "Enzyme_modulation_interaction"

        # --- Multi mechanism  
        if (r1 in PD_REL and r2 in ENZYME_REG) or \
           (r2 in PD_REL and r1 in ENZYME_REG) or \
           ("METABOLIZED_BY" in pair and (r1 in PD_REL or r2 in PD_REL)):
            return "Multi_mechanism"

        return "Target_unclassified"

    elif node == "ATC":
        return "ATC_class"

    elif node == "Metabolite_direct":
        return "PK_metabolite"

    return "Other"
# %%
df_fp["final_mechanism"] = df_fp.apply(derive_final_mechanism, axis=1)
# %%
df_fp = df_fp.merge(
    compound_lookup,
    left_on="herb_compound_key",    
    right_on="compound_key",
    how="left"
).rename(columns={
    "display_name": "herb_compound_name"
}).drop(columns=["compound_key"])

df_fp = df_fp.merge(
    compound_lookup,
    left_on="drug_compound_key",    
    right_on="compound_key",
    how="left"
).rename(columns={
    "display_name": "drug_compound_name"
}).drop(columns=["compound_key"])
 
# %%
df_target = df_fp[df_fp["node_type"] == "Target"].copy()
df_other = df_fp[df_fp["node_type"] != "Target"].copy()
# %%
df_target = df_target.merge(
    target_lookup,
    left_on="node_id",
    right_on="target_id",
    how="left"
)
# %%
df_fp = pd.concat([df_target, df_other], ignore_index=True)
# %%
df_fp.head(100)
# %%
df_fp["display_evidence"] = None
#%%
# Target → ใช้ชื่อจริง
df_fp.loc[df_fp["node_type"] == "Target", "display_evidence"] = df_fp["target_name"].fillna(df_fp["node_id"])

# ATC → ใช้ code
df_fp.loc[df_fp["node_type"] == "ATC", "display_evidence"] = df_fp["node_id"]

# Metabolite → ใส่ label
df_fp.loc[df_fp["node_type"] == "Metabolite_direct", "display_evidence"] = "Direct metabolite link"
#%%
def unique_join(x):
    vals = list(dict.fromkeys(x.dropna()))
    return ", ".join(vals)
summary_df = df_fp.groupby([
    "pair_id", 
    "species_name",  
    "herb_compound_key","drug_compound_key",
    "herb_compound_name","drug_compound_name",
    "probability"
]).agg(
    mechanisms=("final_mechanism", lambda x: ", ".join(sorted(set(x)))),
    shared_evidence=("display_evidence", unique_join)
).reset_index()

final_table = summary_df[[
    "species_name","herb_compound_name","drug_compound_name","probability","mechanisms","shared_evidence"
]]
final_table = final_table.sort_values(
    by=["species_name","drug_compound_name","probability"],
    ascending=[True,True,False]
)
#%%
final_table.head()
#%%
final_table.to_csv("Summary_exFN_table.csv", index=False)
print("Saved: Summary_FN_table.csv")
 # %%



