#%%
import csv
from pathlib import Path
#%%
# -----------------------------
# 1. กำหนด path
# -----------------------------

OUT_DIR = Path("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\drugbank_hdi\\output")

ddi_db_path = OUT_DIR / "ddi_db_pairs.csv"
mapping_path = OUT_DIR / "drugbank_to_chembl.csv"
output_path = OUT_DIR / "ddi_chembl_pairs.csv"

# -----------------------------
# 2. โหลด mapping เป็น dictionary
# -----------------------------

db_to_chembl = {}

with open(mapping_path, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    
    for row in reader:
        db_id = row["drugbank_id"]
        chembl_id = row["chembl_id"]
        
        db_to_chembl[db_id] = chembl_id

print("Total mapping loaded:", len(db_to_chembl))

# -----------------------------
# 3. อ่าน ddi_db_pairs แล้ว map
# -----------------------------

chembl_pairs = set()  # ใช้ set กัน duplicate

with open(ddi_db_path, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    
    for row in reader:
        db1 = row["drugbank_id_1"]
        db2 = row["drugbank_id_2"]
        
        # เช็คว่าทั้งสองตัวมี CHEMBL mapping ไหม
        if db1 in db_to_chembl and db2 in db_to_chembl:
            
            chembl1 = db_to_chembl[db1]
            chembl2 = db_to_chembl[db2]
            
            # ทำ symmetric pair
            pair = tuple(sorted([chembl1, chembl2]))
            
            chembl_pairs.add(pair)

print("Total CHEMBL DDI pairs:", len(chembl_pairs))

# -----------------------------
# 4. Export เป็น CSV
# -----------------------------

with open(output_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    
    writer.writerow(["chembl_id_1", "chembl_id_2"])
    writer.writerows(chembl_pairs)

print("Done.")
# %%
