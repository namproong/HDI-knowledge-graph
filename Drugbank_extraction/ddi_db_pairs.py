#%%
import xml.etree.ElementTree as ET
import csv
from pathlib import Path
#%%
# -------------------------------
# 1. กำหนด path
# -------------------------------

XML_PATH = Path("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\drugbank_hdi\\data\\full database.xml")
OUTPUT_PATH = Path("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\drugbank_hdi\\output\\ddi_db_pairs.csv")

ns = "{http://www.drugbank.ca}"  # namespace ของ DrugBank XML

# set เอาไว้กัน duplicate pair
ddi_pairs = set()

# -------------------------------
# 2. อ่าน XML แบบ streaming
# -------------------------------

for event, elem in ET.iterparse(XML_PATH, events=("end",)):
    
    # --------------------------------
    # 3. เช็คว่า element นี้คือ <drug>
    # --------------------------------
    if elem.tag == f"{ns}drug":
        
        # --------------------------------
        # 4. ดึง primary DrugBank ID ของ drug หลัก
        # --------------------------------
        db_id = elem.findtext(f"{ns}drugbank-id[@primary='true']")
        
        # ถ้าไม่มี primary id ให้ข้าม
        if not db_id:
            elem.clear()
            continue
        
        # --------------------------------
        # 5. เข้าไปดู <drug-interactions>
        # --------------------------------
        for inter in elem.findall(f"{ns}drug-interactions/{ns}drug-interaction"):
            
            # ดึง DrugBank ID ของ drug อีกตัว
            other_id = inter.findtext(f"{ns}drugbank-id")
            
            # ถ้าไม่มี id ให้ข้าม
            if not other_id:
                continue
            
            # --------------------------------
            # 6. ทำให้ pair เป็น symmetric
            # --------------------------------
            pair = tuple(sorted([db_id, other_id]))
            
            # ใส่ลง set (กัน duplicate อัตโนมัติ)
            ddi_pairs.add(pair)
        
        # --------------------------------
        # 7. ล้าง element ออกจาก memory
        # --------------------------------
        elem.clear()

# -------------------------------
# 8. เขียนผลลัพธ์เป็น CSV
# -------------------------------

with open(OUTPUT_PATH, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    
    # เขียน header
    writer.writerow(["drugbank_id_1", "drugbank_id_2"])
    
    # เขียนทุก pair
    writer.writerows(ddi_pairs)

print("Done.")
print("Total DDI pairs:", len(ddi_pairs))
# %%
