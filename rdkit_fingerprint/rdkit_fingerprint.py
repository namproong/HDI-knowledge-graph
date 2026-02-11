#%%
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit import RDLogger
import numpy as np
import pandas as pd
# %%
RDLogger.DisableLog('rdApp.warning')
# %%
morgan_gen = GetMorganGenerator(radius=2, fpSize=2048)
# %%
def mol_to_fp_array(mol):
    fp = morgan_gen.GetFingerprint(mol)
    arr = np.zeros((2048,), dtype=int)
    Chem.DataStructs.ConvertToNumpyArray(fp, arr)
    return arr
# %%
#%%
chembl = pd.read_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data for fp\\chembl_compound_structures_forfinger.csv")

BATCH_SIZE = 2000
out_file = "C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data for fp\\compound_fingerprints_chembl.csv"

first_write = True

for start in range(0, len(chembl), BATCH_SIZE):
    end = start + BATCH_SIZE
    batch = chembl.iloc[start:end]

    rows = []

    for _, row in batch.iterrows():
        if pd.isna(row["molfile"]):
            continue

        mol = Chem.MolFromMolBlock(row["molfile"])
        if mol is None:
            continue

        fp_arr = mol_to_fp_array(mol)

        rows.append({
            "source": "chembl",
            "molregno": row["molregno"],
            "standard_inchi": row["standard_inchi"],
            "standard_inchi_key": row["standard_inchi_key"],
            "fingerprint": fp_arr.tolist()
        })

    if rows:
        pd.DataFrame(rows).to_csv(
            out_file,
            mode="w" if first_write else "a",
            header=first_write,
            index=False
        )
        first_write = False

    print(f"Processed {end}/{len(chembl)}")

# %%
from scipy import sparse
# %%
def inchi_to_morgan(inchi, radius=2, n_bits=2048):
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, radius=radius, nBits=n_bits
        )
        arr = np.zeros((n_bits,), dtype=np.int8)
        Chem.DataStructs.ConvertToNumpyArray(fp, arr)
        return arr
    except:
        return None
# %%
np_df = pd.read_csv(
    "C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data for fp\\natural_products_chembl_merged.csv"
)

BATCH_SIZE = 2000
out_file_np = "C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data for fp\\compound_fingerprints_np.csv"

first_write = True

for start in range(0, len(np_df), BATCH_SIZE):
    end = start + BATCH_SIZE
    batch = np_df.iloc[start:end]

    rows = []

    for _, row in batch.iterrows():
        # เลือก inchi ที่ดีที่สุด
        inchi = row["standard_inchi"]
        if pd.isna(inchi) or inchi in ["n.a.", "NA", ""]:
            inchi = row["InChI"]

        if pd.isna(inchi) or inchi in ["n.a.", "NA", ""]:
            continue

        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            continue

        fp_arr = mol_to_fp_array(mol)

        rows.append({
            "source": "np",
            "np_id": row["np_id"],
            "chembl_id": row.get("chembl_id", None),
            "np_inchi": inchi,
            "np_inchi_key": row.get("standard_inchi_key", None),
            "fingerprint": fp_arr.tolist()
        })

    if rows:
        pd.DataFrame(rows).to_csv(
            out_file_np,
            mode="w" if first_write else "a",
            header=first_write,
            index=False
        )
        first_write = False

    print(f"NP processed {end}/{len(np_df)}")
# %%
def fix_np_inchikey(row):
    if pd.notna(row["np_inchi_key"]):
        return row["np_inchi_key"]
    try:
        return Chem.InchiToInchiKey(row["np_inchi"])
    except:
        return None

np_fp = pd.read_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data for fp\\compound_fingerprints_np.csv")
np_fp["np_inchi_key"] = np_fp.apply(fix_np_inchikey, axis=1)

np_fp.to_csv("C:\\Users\\proon\\OneDrive\\เดสก์ท็อป\\HDi project 222\\data for fp\\compound_fingerprints_np_fixed.csv", index=False)