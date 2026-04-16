#%%
import xml.etree.ElementTree as ET
from pathlib import Path
import csv
# %%
xml_path = Path(
    r"C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data from DB\\full database.xml"
)
# %%
tree = ET.parse(xml_path)
root = tree.getroot()
# %%
ns = {"db": "http://www.drugbank.ca"}
# %%
out_path = Path("drugbank_drug_enzyme_action.csv")
# %%
with open(out_path, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow([
        "drugbank_id",
        "drug_name",
        "chembl_id",
        "enzyme_name",
        "uniprot_id",
        "gene_name",
        "action",
        "organism",
        "tax_id",
        "source"
    ])

    for drug in root.findall("db:drug", ns):

        drugbank_id = drug.findtext(
            "db:drugbank-id[@primary='true']", namespaces=ns
        )
        drug_name = drug.findtext("db:name", namespaces=ns)

        # ChEMBL ID
        chembl_id = None
        for ext in drug.findall(
            "db:external-identifiers/db:external-identifier", ns
        ):
            if ext.findtext("db:resource", namespaces=ns) == "ChEMBL":
                chembl_id = ext.findtext("db:identifier", namespaces=ns)
                break

        # enzymes
        for enz in drug.findall("db:enzymes/db:enzyme", ns):

            enzyme_name = enz.findtext("db:name", namespaces=ns)
            organism = enz.findtext("db:organism", namespaces=ns)

            org_elem = enz.find("db:organism", ns)
            tax_id = (
                org_elem.attrib.get("ncbi-taxonomy-id")
                if org_elem is not None
                else None
            )

            actions = enz.findall("db:actions/db:action", ns)
            action_list = [a.text for a in actions] if actions else ["unknown"]

            poly = enz.find("db:polypeptide", ns)
            if poly is None:
                continue

            uniprot_id = poly.attrib.get("id")
            gene_name = poly.findtext("db:gene-name", namespaces=ns)

            for action in action_list:
                writer.writerow([
                    drugbank_id,
                    drug_name,
                    chembl_id,
                    enzyme_name,
                    uniprot_id,
                    gene_name,
                    action,
                    organism,
                    tax_id,
                    "drugbank"
                ])

print("✅ Exported:", out_path)
# %%
