#%%
import pandas as pd
import numpy as np
import pickle
#%%
#Load Data
raw_data = pd.read_parquet("rotate3con_embeddings_with_fingerprint_only.parquet")
raw_data = raw_data.rename(columns={"node_name": "compound_key"})
print(raw_data.columns)
print(raw_data.shape) #
#%%
# โหลด model bundle
with open("model_bundle.pkl", "rb") as f:
    bundle = pickle.load(f)

# แยกค่าที่จำเป็นออกมา
model = bundle["model"]
threshold = bundle["threshold"]
feature_shape = bundle["feature_shape"]

print("Model loaded!")
print("Threshold:", threshold) #Threshold: 0.40552763819095483
print("Expected number of features:", feature_shape) #1205
#%%
external_pairs = pd.read_csv(
    "drug herb 30 + 20interaction.csv",
    header=None
)
#%%
external_pairs.head()
#%%
external_pairs.columns = [
    "species_tax_id",
    "species_name",
    "compound_key_1",
    "compound_key_2"
]
print(external_pairs.head())
print(external_pairs.shape) #(11925, 4)
# %%
valid_compounds = set(raw_data["compound_key"])
external_pairs = external_pairs[
    external_pairs["compound_key_1"].isin(valid_compounds) &
    external_pairs["compound_key_2"].isin(valid_compounds)
    ].copy()
print("Pairs with embeddings:", external_pairs.shape)
#Pairs with embeddings: (11925, 4)
# %%
embedding_cols = [
    col for col in raw_data.columns 
    if col not in ["compound_key", "fingerprint"]
]
#%%
print(raw_data.columns[:20])
print(len(embedding_cols))
#%%
# Merge Embedding AFTER split
def merge_embedding(pair_df):

    merged = pair_df.merge(
        raw_data,
        left_on="compound_key_1",
        right_on="compound_key",
        how="left"
    ).rename(columns={
        **{col: f"A_{col}" for col in embedding_cols},
        "fingerprint": "fingerprint_A"
    })

    merged = merged.merge(
        raw_data,
        left_on="compound_key_2",
        right_on="compound_key",
        how="left",
        suffixes=("", "_B")
    ).rename(columns={
        **{col: f"B_{col}" for col in embedding_cols},
        "fingerprint": "fingerprint_B"
    })

    merged = merged.drop(columns=["compound_key", "compound_key_B"])

    return merged
#%%
external_merged = merge_embedding(external_pairs)
print(external_merged.shape) #(11925, 806) >> 400 + 400 + cp1 + cp1 + species tax is + species name + fingerprint 1 + fingerprint 1
print(external_merged.filter(regex="A_").shape) #(11925, 400)
print(external_merged.filter(regex="B_").shape) #(11925, 400)
print(external_merged.columns)
#%%
# Feature Engineering 
def compute_tanimoto_chunked(df, chunk_size=100000):

    n = len(df)

    tanimoto = np.zeros(n, dtype=np.float32)
    intersection_arr = np.zeros(n, dtype=np.float32)
    union_arr = np.zeros(n, dtype=np.float32)

    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)

        fp_A_chunk = df["fingerprint_A"].iloc[start:end]
        fp_B_chunk = df["fingerprint_B"].iloc[start:end]

        for i, (fp1, fp2) in enumerate(zip(fp_A_chunk, fp_B_chunk)):

            a = np.fromstring(fp1.strip()[1:-1], sep=',', dtype=np.uint8)
            b = np.fromstring(fp2.strip()[1:-1], sep=',', dtype=np.uint8)

            intersection = np.sum(a & b)
            union = np.sum(a) + np.sum(b) - intersection

            idx = start + i

            intersection_arr[idx] = intersection
            union_arr[idx] = union
            tanimoto[idx] = intersection / union if union != 0 else 0

        print(f"Processed {end} / {n}")

    return np.vstack([
        intersection_arr,
        union_arr,
        tanimoto
    ]).T

def build_features(merged):

    A_cols = [f"A_{col}" for col in embedding_cols]
    B_cols = [f"B_{col}" for col in embedding_cols]

    A = merged[A_cols].values.astype("float32")
    B = merged[B_cols].values.astype("float32")
    hadamard = A * B
    l1 = np.abs(A - B)
    l2 = (A - B) ** 2
    dot = np.sum(A * B, axis=1).reshape(-1,1)
    cos = np.sum(A * B, axis=1) / (
        (np.linalg.norm(A, axis=1) *
        np.linalg.norm(B, axis=1)) + 1e-8
    )
    cos = cos.reshape(-1,1)
    tanimoto_feats = compute_tanimoto_chunked(merged)

    X = np.hstack([
        hadamard,
        l1,
        l2,
        dot,
        cos,
        tanimoto_feats
    ]).astype("float32")

    return X
#%%
X_external = build_features(external_merged)
#%%
print(X_external.shape) #(11925, 1205)
print(model.n_features_in_) #1205
external_prob = model.predict_proba(X_external)[:,1]
external_merged["probability"] = external_prob
external_merged["prediction"] = (external_prob >= threshold).astype(int)
# %%
external_merged.head()
# %%
lookup = pd.read_csv("compound_lookup_master.csv")
# %%
lookup = lookup[["compound_key","display_name"]]
# %%
external_merged = external_merged.merge(
    lookup,
    left_on="compound_key_1",
    right_on="compound_key",
    how="left"
).rename(columns={"display_name":"compound_name_1"}).drop(columns=["compound_key"])

external_merged = external_merged.merge(
    lookup,
    left_on="compound_key_2",
    right_on="compound_key",
    how="left"
).rename(columns={"display_name":"compound_name_2"}).drop(columns=["compound_key"])
# %%
# Build Compound-Level Evidence Table
compound_evidence = external_merged[
[
    "species_tax_id",
    "species_name",
    "compound_key_1",
    "compound_name_1",
    "compound_key_2",
    "compound_name_2",
    "probability",
    "prediction"
]]
#%%
compound_evidence.head(50)
# %%
compound_evidence.to_csv(
    "compound_level_evidence.csv",
    index=False
)
# %%
predicted_pairs = compound_evidence[
    compound_evidence["prediction"] == 1
].copy()
# %%
# aggregate จาก compound → herb-drug level
herb_drug_summary = (
    predicted_pairs
    .groupby(["species_name", "compound_name_2"])
    .agg(
        supporting_compounds=("compound_name_1",
            lambda x: "; ".join(sorted(set(x)))),
        max_probability=("probability", "max")
    )
    .reset_index()
)

herb_drug_summary = herb_drug_summary.rename(columns={
    "species_name": "herb_species",
    "compound_name_2": "drug_name"
})

herb_drug_summary = herb_drug_summary.sort_values(
    ["max_probability"],
    ascending=[False]
)

print(herb_drug_summary.head(50))
# %%
len(herb_drug_summary) #48
# %%
#รวม compound ใน herb 
external_herb_drug = (
    compound_evidence
    .groupby(["species_name", "compound_name_2"])
    .size()
    .reset_index(name="reported_compounds") #แปลว่าไร ใส่ทำไม
)

external_herb_drug = external_herb_drug.rename(columns={
    "species_name": "herb_species",
    "compound_name_2": "drug_name"
})

print("External reported herb-drug pairs:", external_herb_drug.shape) #(54, 3)
external_herb_drug.head(100) 
# %%
#%%Compare reported vs predicted
evaluation = external_herb_drug.merge(
    herb_drug_summary,
    on=["herb_species","drug_name"],
    how="left"
)
#%%
# เรียงตาม confidence
evaluation = evaluation.sort_values(
    ["max_probability"],
    ascending=[False]
)
evaluation.head(100)
# %%
len(evaluation) #54
# %%
#%% 
evaluation_summary = evaluation[[
    "herb_species",
    "drug_name",
    "max_probability",
    "supporting_compounds"

]]


evaluation_summary.to_csv(
    "all_external_hdi_evaluation.csv",
    index=False
)
# %%

# coverage analysis
total_pairs = len(evaluation)
supported_pairs = len(herb_drug_summary) 
unsupported_pairs = total_pairs - supported_pairs

coverage = supported_pairs / total_pairs

print("Total reported herb-drug pairs:", total_pairs) #54
print("Supported by model:", supported_pairs) #48
print("Unsupported:", unsupported_pairs) #6
print("Coverage:", round(coverage*100,2), "%") #88.89 %
# %%
tp_compound_evidence = compound_evidence[
    compound_evidence["prediction"] == 1
].copy()

# เรียงตาม species
tp_compound_evidence = tp_compound_evidence.sort_values(
    ["species_name", "compound_name_1", "compound_name_2"]
)

tp_compound_evidence.to_csv("exTP_compound_level_evidence.csv", index=False)
# %%
tp_compound_evidence.head(100)
# %%
num_unique_pairs = tp_compound_evidence[
    ["species_name", "compound_name_2"]
].drop_duplicates().shape[0]

print(num_unique_pairs) #48 
# %%
# รายชื่อ species-drug ที่ FN
fn_pairs = [
    ("Cymbopogon citratus", "CANFOSFAMIDE"),
    ("Hypericum perforatum", "CYCLOSPORINE"),
    ("Lessertia frutescens", "ATAZANAVIR SULFATE"),
    ("Peumus boldus", "WARFARIN POTASSIUM"),
    ("Piper nigrum", "THEOPHYLLINE"),
    ("Salix fragilis", "EDOXABAN")
]

fn_compound_evidence = pd.DataFrame()

for herb, drug in fn_pairs:
    df = compound_evidence[
        (compound_evidence["species_name"] == herb) &
        (compound_evidence["compound_name_2"] == drug)
    ].copy()
    fn_compound_evidence = pd.concat([fn_compound_evidence, df], axis=0)

fn_compound_evidence = fn_compound_evidence.sort_values(
    ["species_name", "compound_name_1", "compound_name_2"]
)

fn_compound_evidence.to_csv("exFN_compound_level_evidence.csv", index=False)
# %%
total_rows = len(fn_compound_evidence)
print("Total FN rows:", total_rows) 
num_unique_pairs = fn_compound_evidence[
    ["species_name", "compound_name_2"]
].drop_duplicates().shape[0]

print("Unique species-drug pairs (FN):", num_unique_pairs)
rows_per_species = fn_compound_evidence.groupby("species_name").size()

print(rows_per_species)
# %%
 