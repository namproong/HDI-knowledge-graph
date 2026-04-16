#%%
import pandas as pd
# %%
# load data 
species_compound = pd.read_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\species_compound_notF.csv")
compound_lookup = pd.read_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\compound_lookup_master.csv")
plant_compound = pd.read_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\plant_compound.csv")
species_info_selected = pd.read_csv("C:\\Users\\proon\\Downloads\\species_info_selected.csv")
#%%
# merge (inner join ตัด compound ที่ไม่มี display_name ออกเพราะเป็น compound ที่ไม่อยู่ในระบบ)
species_compound_display = species_compound.merge(
    compound_lookup[["compound_key", "display_name"]],
    on="compound_key",
    how="inner"
)

# เรียงเพื่อให้อ่านง่าย
species_compound_display = species_compound_display.sort_values(
    ["species_tax_id", "compound_key"]
)
# %%
species_compound_display.columns #Index(['species_tax_id', 'species_name', 'compound_key', 'display_name'], dtype='object')
# %%
species_compound_display.head()
# %%
species_compound_display["species_tax_id"].count() #644118 
# %%
len(species_compound_display) #644118 
# %%
species_compound_display["compound_key"].count()  #644118  # %%
# %%
#%%
plant_species = species_info_selected[
    species_info_selected["kingdom_name"] == "Viridiplantae"
]
# %%
plant_species_ids = plant_species["species_tax_id"].drop_duplicates()
# %%
print(species_compound_display["species_tax_id"].dtype)
print(plant_species_ids.dtype)
#%%
species_compound_display["species_tax_id"] = species_compound_display["species_tax_id"].astype(str)
plant_species_ids = plant_species_ids.astype(str)
#%%
species_compound_display = species_compound_display[
    species_compound_display["species_tax_id"].isin(plant_species_ids)
] #เก็บเฉพาะแถวที่ species_tax_id อยู่ใน list ของ plant species
# %%
print("rows:", len(species_compound_display)) #475359
print("species:", species_compound_display["species_tax_id"].nunique()) #14206 
print("compound:", species_compound_display["compound_key"].nunique()) #100383
# %%
species_compound_display.to_csv(
    "C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\species_compound_display.csv",
    index=False
)
# %%
