#%%
import pandas as pd
import random
import numpy as np
#%% 
from sklearn.model_selection import train_test_split
#%% 
random.seed(42)
np.random.seed(42)
#%%
#Load Data
raw_data = pd.read_parquet("rotate3con_embeddings_with_fingerprint_only.parquet")
raw_data = raw_data.rename(columns={"node_name": "compound_key"})
print(raw_data.columns)
#%%
ddi = pd.read_csv("HDI_ddi_compound_pairs.csv")
#%%
print(raw_data.shape) #208064, 402
print(ddi.shape) #1059416, 2
#%% 
# Embedding universe 
valid_compounds = set(raw_data["compound_key"])
print(len(valid_compounds)) #208064
#%% 
# จำกัด pairs ให้อยู่ใน embedding universe
pairs_filtered = ddi[
    ddi["compound_key_1"].isin(valid_compounds) &
    ddi["compound_key_2"].isin(valid_compounds)
].copy()
print(pairs_filtered.shape) #(912547, 2)
#%%
# Normalize pair กัน order leakage 
def normalize_pair(a, b):
    return tuple(sorted([a, b]))

pairs_filtered["pair"] = [
    normalize_pair(a, b)
    for a, b in zip(
        pairs_filtered["compound_key_1"],
        pairs_filtered["compound_key_2"]
    )
]

pairs_filtered = pairs_filtered.drop_duplicates("pair")
print(pairs_filtered.shape) #(912547, 3)
#%%
print(normalize_pair("A", "B"))
print(normalize_pair("B", "A"))
#%%
# Train / Test Split (Pair-level split)
train_pairs, test_pairs = train_test_split(
    pairs_filtered["pair"].values,
    test_size=0.2,
    random_state=42
)

train_pos_pairs = set(train_pairs)
test_pos_pairs = set(test_pairs)

print("Train positives:", len(train_pos_pairs)) #Train positives: 730037
print("Test positives:", len(test_pos_pairs)) #Test positives: 182510
#%%
set(train_pos_pairs).intersection(set(test_pos_pairs)) #0
#%%
# Generate Negative Pairs (separate train/test)
all_positive = train_pos_pairs.union(test_pos_pairs)

graph_compounds = set(
    pairs_filtered["pair"].apply(lambda x: x[0])
).union(
    pairs_filtered["pair"].apply(lambda x: x[1])
)

all_compounds = list(graph_compounds)

def generate_negative_pairs(n_samples, forbidden_pairs):
    negative_pairs = set()

    while len(negative_pairs) < n_samples:
        a, b = random.sample(all_compounds, 2)
        pair = normalize_pair(a, b)

        if pair not in forbidden_pairs and pair not in negative_pairs:
            negative_pairs.add(pair)

    return negative_pairs

train_neg_pairs = generate_negative_pairs(
    len(train_pos_pairs),
    all_positive
)
test_forbidden = all_positive.union(train_neg_pairs)

test_neg_pairs = generate_negative_pairs(
    len(test_pos_pairs),
    test_forbidden
)

print("Train Negatives:", len(train_neg_pairs)) #Train Negatives: 730037
print("Test Negatives:", len(test_neg_pairs)) #Test Negatives: 182510
print(
    "Train/Test negative overlap:",
    len(train_neg_pairs & test_neg_pairs)) #Train/Test negative overlap: 0 
#%%
# Build DataFrame
def build_df(pos_pairs, neg_pairs):
    pos_df = pd.DataFrame(list(pos_pairs), columns=["compound_key_1","compound_key_2"])
    pos_df["label"] = 1

    neg_df = pd.DataFrame(list(neg_pairs), columns=["compound_key_1","compound_key_2"])
    neg_df["label"] = 0

    df_all = pd.concat([pos_df, neg_df], ignore_index=True)
    df_all = df_all.sample(frac=1, random_state=42).reset_index(drop=True)

    return df_all

train_final = build_df(train_pos_pairs, train_neg_pairs)
test_final = build_df(test_pos_pairs, test_neg_pairs)
#%%
set(map(tuple, train_final[["compound_key_1","compound_key_2"]].values)) & set(map(tuple, test_final[["compound_key_1","compound_key_2"]].values))
print(len(train_final), len(test_final)) #(1460074, 365020)
print(train_final.shape) #(1460074, 3)
print(test_final.shape) #(365020, 3) 
# %%
train_final.to_csv("train_pairsHDIfinal.csv", index=False)
test_final.to_csv("test_pairsHDIfinal.csv", index=False)
#%%
print(raw_data.dtypes.head())
#%%
embedding_cols = [
    col for col in raw_data.columns 
    if col not in ["compound_key", "fingerprint"]
]

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
train_merged = merge_embedding(train_final)
#%%
test_merged = merge_embedding(test_final)
#%%
print(train_merged.shape) #(1460074, 805)
print(test_merged.shape) #(365020, 805)
print(train_merged.filter(regex="A_").shape) #(1460074, 400)
print(train_merged.filter(regex="B_").shape) #(1460074, 400)
#%%
print(train_merged.columns[:50])
print(test_merged.columns[:50])
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
    y = merged["label"].values

    return X, y
#%%
X_train, y_train = build_features(train_merged)
#%%
X_test, y_test = build_features(test_merged)
#%%
print(X_train.shape) #(1460074, 1205)
print(X_test.shape) # (365020, 1205)
#%%
np.save("X_trainHDIfinal.npy", X_train.astype("float32"))
np.save("y_trainHDIfinal.npy", y_train.astype("int8"))
#%%
np.save("X_testHDIfinal.npy", X_test.astype("float32"))
np.save("y_testHDIfinal.npy", y_test.astype("int8"))
#%%
