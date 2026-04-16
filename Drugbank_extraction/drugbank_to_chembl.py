#%%
import xml.etree.ElementTree as ET
import csv
#%%
from pathlib import Path
#%%
XML_PATH = Path("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\drugbank_hdi\\data\\full database.xml")
OUT_DIR = Path("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\drugbank_hdi\\output")
OUT_DIR.mkdir(exist_ok=True)
OUTPUT_PATH = OUT_DIR / "drugbank_to_chembl.csv"

ns = "{http://www.drugbank.ca}"

chembl_map = []

for event, elem in ET.iterparse(XML_PATH, events=("end",)):
    
    if elem.tag == f"{ns}drug":
        
        # --- 1. get primary DrugBank ID ---
        db_id = elem.findtext(f"{ns}drugbank-id[@primary='true']")
        
        # --- 2. find ChEMBL id ---
        chembl_id = None
        
        for ext in elem.findall(f"{ns}external-identifiers/{ns}external-identifier"):
            
            resource = ext.findtext(f"{ns}resource")
            identifier = ext.findtext(f"{ns}identifier")
            
            if resource == "ChEMBL":
                chembl_id = identifier
                break
        
        # --- 3. save if found ---
        if db_id and chembl_id:
            chembl_map.append((db_id, chembl_id))
        
        elem.clear()

# --- 4. export CSV ---
with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["drugbank_id", "chembl_id"])
    writer.writerows(chembl_map)

print("Done.")
print("Total mapped:", len(chembl_map))
# %%
