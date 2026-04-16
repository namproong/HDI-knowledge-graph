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
# %%
# load prediction results
test_pairs = pd.read_csv("final_test_out.csv")
print("Test pairs:", test_pairs.shape) #Test pairs: (365020, 5) 
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
# %%
TP = test_pairs[(test_pairs.y_true==1) & (test_pairs.prediction==1)]
FP = test_pairs[(test_pairs.y_true==0) & (test_pairs.prediction==1)]
FN = test_pairs[(test_pairs.y_true==1) & (test_pairs.prediction==0)]
TN = test_pairs[(test_pairs.y_true==0) & (test_pairs.prediction==0)]

print("TP:", len(TP)) #TP: 156891
print("FP:", len(FP)) #FP: 66409
print("FN:", len(FN)) #FN: 25619
print("TN:", len(TN)) #TN: 116101
#%%
top_tp = TP.sort_values("probability", ascending=False).head(100)
print(top_tp.head()) #compound1 compound2 y_true  probability  prediction 
top_tp.to_csv("final_top_true_positive_pairs.csv", index=False)
print("Saved: top_true_positive_pairs.csv")
# %%
#list compound pairs สำหรับ Neo4j
tp_pairs = list(zip(top_tp.compound_key_1, top_tp.compound_key_2))
print("TP pairs:", len(tp_pairs)) #TP pairs: 100
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

tp_evidence_rows = []

with driver.session(database="hdikgv3") as session:
    for i, row in top_tp.iterrows():

        pair_id = i
        c1 = row["compound_key_1"]
        c2 = row["compound_key_2"]
        prob = row["probability"]

        for query in all_queries:

            results = session.execute_read(run_query, query, c1, c2)

            for row_q in results:
                tp_evidence_rows.append([
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

for row in tp_evidence_rows:

    pair_id = row[0]
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
    "compound_key_1",
    "compound_key_2",
    "probability",
    "node_type",
    "node_id",
    "rel_c1",
    "rel_c2",
    "mechanism_c1",
    "mechanism_c2"
]
df = pd.DataFrame(final_rows, columns=columns)
#%%
df["node_id"] = df["node_id"].astype(str)
target_lookup["target_id"] = target_lookup["target_id"].astype(str)
#%%
priority = {
    "Target": 1,
    "Metabolite_direct": 2,
    "ATC": 3
}
df["priority"] = df["node_type"].map(priority)
# sort ตาม pair + ความสำคัญ
df = df.sort_values(
    by=["pair_id", "priority"],
)
#%%
df = df.drop(columns=["priority"])
df = df.drop_duplicates(
    subset=[
        "pair_id",
        "compound_key_1",
        "compound_key_2",
        "node_type",
        "node_id",
        "rel_c1",
        "rel_c2"
    ]
)
#%%
df.head(200)
#%%
df.to_csv("final_TP_evidence_table.csv", index=False)
print("Saved: final_TP_evidence_table.csv")
#%%
df.head(100)
#%%
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
#%%
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
#%%
df["final_mechanism"] = df.apply(derive_final_mechanism, axis=1)
#%%
df = df.merge(
    compound_lookup,
    left_on="compound_key_1",
    right_on="compound_key",
    how="left"
).rename(columns={
    "display_name": "compound_name_1"
}).drop(columns=["compound_key"])
#%%
df = df.merge(
    compound_lookup,
    left_on="compound_key_2",
    right_on="compound_key",
    how="left"
).rename(columns={
    "display_name": "compound_name_2"
}).drop(columns=["compound_key"])

#%%
df_target = df[df["node_type"] == "Target"].copy()
df_other = df[df["node_type"] != "Target"].copy()
#%%
df_target = df_target.merge(
    target_lookup,
    left_on="node_id",
    right_on="target_id",
    how="left"
)
#%%
df = pd.concat([df_target, df_other], ignore_index=True)
#%%
df.head(100)
#%%
df["display_evidence"] = None
#%%
# Target → ใช้ชื่อจริง
df.loc[df["node_type"] == "Target", "display_evidence"] = df["target_name"].fillna(df["node_id"])

# ATC → ใช้ code
df.loc[df["node_type"] == "ATC", "display_evidence"] = df["node_id"]

# Metabolite → ใส่ label
df.loc[df["node_type"] == "Metabolite_direct", "display_evidence"] = "Direct metabolite link"
#%%
def unique_join(x):
    vals = list(dict.fromkeys(x.dropna()))
    return ", ".join(vals)

summary_df = (
    df
    .groupby([
        "pair_id",
        "compound_key_1",
        "compound_key_2",
        "compound_name_1",
        "compound_name_2",
        "probability"
    ])
    .agg(
        mechanisms=("final_mechanism",
            lambda x: ", ".join(sorted(set(x)))),

        shared_evidence=("display_evidence",
            lambda x: unique_join(x))
    )
    .reset_index()
)
#%%
final_table = summary_df[[
    "compound_name_1",
    "compound_name_2",
    "probability",
    "mechanisms",
    "shared_evidence"
]]
#%%
final_table = final_table.sort_values(
    by="probability",
    ascending=False
)
final_table.to_csv("Summary_TP_evidence_table.csv", index=False)
print("Saved: Summary_TP_table.csv")
# %%
df[df["final_mechanism"] == "Other"][[
    "rel_c1", "rel_c2"
]].drop_duplicates()
# %%
df[df["final_mechanism"] == "Target_unclassified"][[
    "rel_c1", "rel_c2"
]].drop_duplicates()
#%%
df[df["node_type"] == "Target"][[
    "node_id", "target_name"
]].drop_duplicates().head(20)
# %%
evidence_count = df.groupby([
    "compound_key_1",
    "compound_key_2"
]).size().reset_index(name="n_evidence")
# %%
top_tp_pairs = top_tp[[
    "compound_key_1",
    "compound_key_2"
]]

check = top_tp_pairs.merge(
    evidence_count,
    on=["compound_key_1", "compound_key_2"],
    how="left"
)
#%%
no_evidence = check[check["n_evidence"].isna()]
print("Pairs with NO evidence:", len(no_evidence))
# %%
