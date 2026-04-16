#%% 
import pandas as pd
#%% 
# load model dataset
df_model = pd.read_parquet("C:\\Users\\proon\\Downloads\\final_data_for_xgboost.parquet")
# %%
# load display table
df_display = pd.read_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\n_unified_display_name.csv")
# %%
print(df_model.shape) #(310694, 411)
print(df_display.shape) # (225019, 3)
# %%
df_model = df_model[df_model["fingerprint"].notnull()]
# %%
print(df_model.shape) #(208064, 411)
# %%
df_merged = df_model.merge(
    df_display,
    on="compound_key",
    how="left"
)

print(df_merged.shape) (208064, 413) 
print("Missing display:", df_merged["display_name"].isna().sum()) #Missing display: 1082
# %%
df_merged.columns
# %%
missing_keys = df_merged.loc[
    df_merged["display_name"].isna(),
    "compound_key"
]

print(missing_keys.head())
print(len(missing_keys))
# %%
df_merged["display_name"] = df_merged["display_name"].fillna(
    df_merged["compound_key"]
)
# %%
df_merged.columns
# %%
df_merged[
    df_merged["compound_key"] == "CMP_IIQKUGXEGMZCLE-UHFFFAOYSA-N"
]
# %%
mask = df_merged["display_name"] == df_merged["compound_key"]

df_merged.loc[mask, "display_source"] = "compound_key_fallback"
# %%
print("Missing display_name:", df_merged["display_name"].isna().sum())
print("Missing display_source:", df_merged["display_source"].isna().sum())
# %%
df_merged[
    df_merged["compound_key"] == "CMP_IIQKUGXEGMZCLE-UHFFFAOYSA-N"
]
# %%
compound_lookup = df_merged[
    ["compound_key", "display_name", "display_source"]
].drop_duplicates()

print(compound_lookup.shape) #(208064, 3)
# %%
compound_lookup["compound_key"].nunique()
# %%
compound_lookup.to_csv(
    "C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\compound_lookup_master.csv",
    index=False
)
# %%
df_model.columns
# %%
