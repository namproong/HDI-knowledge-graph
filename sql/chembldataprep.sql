CREATE SCHEMA IF NOT EXISTS hdi;
SET search_path TO hdi;

CREATE TABLE hdi_stg_chembl_fp (
    source                TEXT,
    molregno              INTEGER,
    standard_inchi        TEXT,
    standard_inchi_key    VARCHAR(27),
    fingerprint           JSONB
);

CREATE TABLE hdi_stg_np_fp (
    source                TEXT,
    np_id                 VARCHAR(100),
    chembl_id             VARCHAR(50),
    np_inchi	          TEXT,
    np_inchi_key          VARCHAR(27),
    fingerprint           JSONB
);

CREATE TABLE hdi.hdi_compound (
    compound_key       VARCHAR(40) PRIMARY KEY,
    compound_inchi      TEXT,
    compound_inchi_key  VARCHAR(27) UNIQUE
);

INSERT INTO hdi.hdi_compound (
    compound_key,
    compound_inchi,
    compound_inchi_key
)
SELECT DISTINCT
    'CMP_' || inchi_key AS compound_key,
    inchi,
    inchi_key
FROM (
    -- ChEMBL (standard)
    SELECT
        standard_inchi      AS inchi,
        standard_inchi_key  AS inchi_key
    FROM hdi.hdi_stg_chembl_fp
    WHERE standard_inchi_key IS NOT NULL

    UNION

    -- Natural Product (used / fallback)
    SELECT
        np_inchi           AS inchi,
        np_inchi_key       AS inchi_key
    FROM hdi.hdi_stg_np_fp
    WHERE np_inchi_key IS NOT NULL
) t;

ALTER TABLE hdi.hdi_stg_chembl_fp
ADD COLUMN chembl_id VARCHAR(50);

UPDATE hdi.hdi_stg_chembl_fp fp
SET chembl_id = md.chembl_id
FROM public.molecule_dictionary md
WHERE fp.molregno = md.molregno;

CREATE TABLE hdi.hdi_compound_mapping (
    compound_key   VARCHAR(40) REFERENCES hdi.hdi_compound(compound_key),
    source         TEXT,                 -- 'chembl', 'np'
    source_id      TEXT,                 -- chembl_id / np_id
    source_sub_id  TEXT,                 -- molregno (ถ้ามี)
    PRIMARY KEY (compound_key, source, source_id)
);

INSERT INTO hdi.hdi_compound_mapping (
    compound_key,
    source,
    source_id,
    source_sub_id
)
SELECT DISTINCT
    c.compound_key,
    'chembl' AS source,
    fp.chembl_id AS source_id,
    fp.molregno::TEXT AS source_sub_id
FROM hdi.hdi_compound c
JOIN hdi.hdi_stg_chembl_fp fp
  ON c.compound_inchi_key = fp.standard_inchi_key;


INSERT INTO hdi.hdi_compound_mapping
SELECT DISTINCT
    c.compound_key,
    'np' AS source,
    fp.np_id AS source_id,
    fp.chembl_id AS source_sub_id
FROM hdi.hdi_compound c
JOIN hdi.hdi_stg_np_fp fp
  ON c.compound_inchi_key = fp.np_inchi_key;

CREATE TABLE hdi.hdi_compound_properties (
    compound_key VARCHAR(40) PRIMARY KEY
        REFERENCES hdi.hdi_compound(compound_key),

    mw_freebase NUMERIC,
    full_mwt NUMERIC,
    alogp NUMERIC,
    hba INTEGER,
    hbd INTEGER,
    psa NUMERIC,
    rtb INTEGER,
    aromatic_rings INTEGER,
    heavy_atoms INTEGER,
    num_ro5_violations INTEGER,
    qed_weighted NUMERIC,
    np_likeness_score NUMERIC
);
INSERT INTO hdi.hdi_compound_properties (
    compound_key,
    mw_freebase,
    full_mwt,
    alogp,
    hba,
    hbd,
    psa,
    rtb,
    aromatic_rings,
    heavy_atoms,
    num_ro5_violations,
    qed_weighted,
    np_likeness_score
)
SELECT DISTINCT ON (c.compound_key)
    c.compound_key,
    cp.mw_freebase,
    cp.full_mwt,
    cp.alogp,
    cp.hba,
    cp.hbd,
    cp.psa,
    cp.rtb,
    cp.aromatic_rings,
    cp.heavy_atoms,
    cp.num_ro5_violations,
    cp.qed_weighted,
    cp.np_likeness_score
FROM hdi.hdi_compound c
JOIN hdi.hdi_stg_chembl_fp fp
  ON c.compound_inchi_key = fp.standard_inchi_key
JOIN public.molecule_dictionary md
  ON fp.molregno = md.molregno
LEFT JOIN public.compound_properties cp
  ON md.molregno = cp.molregno
ORDER BY c.compound_key;
INSERT INTO hdi.hdi_compound_properties (compound_key)
SELECT compound_key
FROM hdi.hdi_compound
ON CONFLICT (compound_key) DO NOTHING;

CREATE TABLE hdi.hdi_compound_name (
    compound_key VARCHAR(40)
        REFERENCES hdi.hdi_compound(compound_key),
    name_type TEXT CHECK (name_type IN ('pref_name', 'synonym')),
    name TEXT NOT NULL,
    source TEXT DEFAULT 'chembl',
    PRIMARY KEY (compound_key, name, name_type)
);
INSERT INTO hdi.hdi_compound_name (
    compound_key,
    name_type,
    name,
    source
)
SELECT DISTINCT
    m.compound_key,
    'pref_name',
    md.pref_name,
    'chembl'
FROM hdi.hdi_compound_mapping m
JOIN public.molecule_dictionary md
  ON md.molregno = m.source_sub_id::INTEGER
WHERE m.source = 'chembl'
  AND md.pref_name IS NOT NULL;
 
INSERT INTO hdi.hdi_compound_name (
    compound_key,
    name_type,
    name,
    source
)
SELECT
    m.compound_key,
    'synonym',
    ms.synonyms,
    'chembl'
FROM hdi.hdi_compound_mapping m
JOIN public.molecule_synonyms ms
  ON ms.molregno = m.source_sub_id::INTEGER
WHERE m.source = 'chembl'
  AND ms.synonyms IS NOT NULL
GROUP BY
    m.compound_key,
    ms.synonyms;

CREATE TABLE hdi.hdi_compound_natural_product (
    compound_key VARCHAR(40)
        REFERENCES hdi.hdi_compound(compound_key),
    source TEXT CHECK (source IN ('chembl', 'external_np')),
    PRIMARY KEY (compound_key, source)
);
INSERT INTO hdi.hdi_compound_natural_product (compound_key, source)
SELECT DISTINCT
    m.compound_key,
    'chembl'
FROM hdi.hdi_compound_mapping m
JOIN public.molecule_dictionary md
  ON md.molregno = m.source_sub_id::INTEGER
WHERE m.source = 'chembl'
  AND md.natural_product = 1;

INSERT INTO hdi.hdi_compound_natural_product (compound_key, source)
SELECT DISTINCT
    compound_key,
    'external_np'
FROM hdi.hdi_compound_mapping
WHERE source = 'np';

CREATE TABLE hdi.hdi_compound_fingerprint (
    compound_key VARCHAR(40)
        REFERENCES hdi.hdi_compound(compound_key),
    fp_type TEXT NOT NULL,          
    fingerprint JSONB NOT NULL,    
    fp_source TEXT,                 
    PRIMARY KEY (compound_key, fp_type)
);
INSERT INTO hdi.hdi_compound_fingerprint (
    compound_key,
    fp_type,
    fingerprint,
    fp_source
)
SELECT
    c.compound_key,
    'ECFP4',
    fp.fingerprint,
    'chembl'
FROM hdi.hdi_compound c
JOIN hdi.hdi_stg_chembl_fp fp
  ON c.compound_inchi_key = fp.standard_inchi_key
WHERE fp.fingerprint IS NOT NULL;

INSERT INTO hdi.hdi_compound_fingerprint (
    compound_key,
    fp_type,
    fingerprint,
    fp_source
)
SELECT
    c.compound_key,
    'ECFP4',
    fp.fingerprint,
    'np'
FROM hdi.hdi_compound c
JOIN hdi.hdi_stg_np_fp fp
  ON c.compound_inchi_key = fp.np_inchi_key
WHERE fp.fingerprint IS NOT NULL
ON CONFLICT (compound_key, fp_type) DO NOTHING;

CREATE TABLE hdi.hdi_target (
    target_key       BIGSERIAL PRIMARY KEY,
    source           TEXT NOT NULL,              -- 'chembl'
    source_target_id TEXT NOT NULL,              -- CHEMBL_ID
    tid              INT,                        -- original TID (debug / join)
    pref_name        TEXT NOT NULL,
    target_type      TEXT,
    organism         TEXT,
    tax_id           INT,
    is_protein       BOOLEAN,
    UNIQUE (source, source_target_id)
);

INSERT INTO hdi.hdi_target (
    source,
    source_target_id,
    tid,
    pref_name,
    target_type,
    organism,
    tax_id,
    is_protein
)
SELECT
    'chembl',
    td.chembl_id,
    td.tid,
    td.pref_name,
    td.target_type,
    td.organism,
    td.tax_id,
    CASE
        WHEN tt.parent_type = 'PROTEIN' THEN TRUE
        ELSE FALSE
    END AS is_protein
FROM target_dictionary td
LEFT JOIN target_type tt
  ON td.target_type = tt.target_type;

CREATE TABLE hdi.hdi_target_protein_class (
    target_key       BIGINT REFERENCES hdi.hdi_target,
    protein_class_id INT,
    class_level      INT,
    pref_name        TEXT,
    PRIMARY KEY (target_key, protein_class_id)
);

INSERT INTO hdi.hdi_target_protein_class (
    target_key,
    protein_class_id,
    class_level,
    pref_name
)
SELECT DISTINCT
    ht.target_key,
    pc.protein_class_id,
    pc.class_level,
    pc.pref_name
FROM hdi.hdi_target ht
JOIN target_components tc
  ON ht.tid = tc.tid
JOIN component_class cc
  ON tc.component_id = cc.component_id
JOIN protein_classification pc
  ON cc.protein_class_id = pc.protein_class_id;

CREATE TABLE hdi.hdi_protein_class_hierarchy (
    protein_class_id INT,
    parent_id  INT,
    PRIMARY KEY (protein_class_id, parent_id)
);

INSERT INTO hdi.hdi_protein_class_hierarchy (
    protein_class_id,
    parent_id
)
SELECT
    protein_class_id,
    parent_id
FROM protein_classification
WHERE parent_id IS NOT NULL;

CREATE TABLE hdi.hdi_target_subset (
    child_target_key  BIGINT NOT NULL,
    parent_target_key BIGINT NOT NULL,
    PRIMARY KEY (child_target_key, parent_target_key),

    CONSTRAINT fk_subset_child
        FOREIGN KEY (child_target_key)
        REFERENCES hdi.hdi_target (target_key),

    CONSTRAINT fk_subset_parent
        FOREIGN KEY (parent_target_key)
        REFERENCES hdi.hdi_target (target_key),

    CONSTRAINT chk_not_self_subset
        CHECK (child_target_key <> parent_target_key)
);
INSERT INTO hdi.hdi_target_subset (
    child_target_key,
    parent_target_key
)
SELECT DISTINCT
    CASE
        WHEN tr.relationship = 'SUBSET OF'
            THEN child.target_key
        WHEN tr.relationship = 'SUPERSET OF'
            THEN parent.target_key
    END AS child_target_key,

    CASE
        WHEN tr.relationship = 'SUBSET OF'
            THEN parent.target_key
        WHEN tr.relationship = 'SUPERSET OF'
            THEN child.target_key
    END AS parent_target_key

FROM target_relations tr

-- target ที่เป็น TID หลัก
JOIN hdi.hdi_target child
  ON tr.tid = child.tid

-- target ที่เป็น RELATED_TID
JOIN hdi.hdi_target parent
  ON tr.related_tid = parent.tid

WHERE tr.relationship IN ('SUBSET OF', 'SUPERSET OF')

ON CONFLICT DO NOTHING;

SELECT
    molregno,
    COUNT(DISTINCT record_id) AS n_record
FROM compound_records
GROUP BY molregno
HAVING COUNT(DISTINCT record_id) > 1
ORDER BY n_record DESC
LIMIT 20;
SELECT
    dm.record_id,
    dm.molregno,
    dm.mechanism_of_action,
    dm.action_type,
    td.pref_name AS target_name
FROM drug_mechanism dm
LEFT JOIN target_dictionary td
  ON dm.tid = td.tid
WHERE dm.molregno = 241
ORDER BY dm.record_id;

CREATE TABLE hdi.hdi_compound_target_pd (
    compound_key        VARCHAR(40)
        REFERENCES hdi.hdi_compound(compound_key),
    target_key          BIGINT,
    action_type         TEXT,

    direct_interaction  BOOLEAN,
    molecular_mechanism BOOLEAN,
    disease_efficacy    BOOLEAN,

    source              TEXT,
    evidence_count      INT,
    source_ids          TEXT[],

    PRIMARY KEY (compound_key, target_key, action_type)
);
INSERT INTO hdi.hdi_compound_target_pd (
    compound_key,
    target_key,
    action_type,
    direct_interaction,
    molecular_mechanism,
    disease_efficacy,
    source,
    evidence_count,
    source_ids
)
SELECT
    cm.compound_key,
    t.target_key,
    dm.action_type,

    BOOL_OR(dm.direct_interaction = 1)  AS direct_interaction,
    BOOL_OR(dm.molecular_mechanism = 1) AS molecular_mechanism,
    BOOL_OR(dm.disease_efficacy = 1)    AS disease_efficacy,

    'chembl'                            AS source,
    COUNT(*)                            AS evidence_count,
    ARRAY_AGG(dm.mec_id::TEXT)          AS source_ids
FROM drug_mechanism dm
JOIN compound_records cr
       ON dm.record_id = cr.record_id
JOIN hdi.hdi_compound_mapping cm
       ON cm.source = 'chembl'
      AND cm.source_sub_id = cr.molregno::TEXT
JOIN hdi.hdi_target t
       ON dm.tid = t.tid
WHERE dm.tid IS NOT NULL
GROUP BY
    cm.compound_key,
    t.target_key,
    dm.action_type;

CREATE TABLE hdi.hdi_stg_final_target_chembl_data (
    target_id    TEXT,
    target_name  TEXT,
    target_type  TEXT,
    chembl_id    TEXT
);
SELECT COUNT(*) FROM hdi_stg_final_target_chembl_data;
CREATE TABLE stg_activities_np_target_selected (
    np_id     TEXT,
    target_id TEXT
);

INSERT INTO hdi.hdi_compound_target_pd (
    compound_key,
    target_key,
    source,
    evidence_count,
    source_ids
)
SELECT
    cm.compound_key,
    t.target_key,
    'np' AS source,
    COUNT(*) AS evidence_count,
	NULL::TEXT[] AS source_ids
FROM hdi.stg_activities_np_target_selected s
JOIN hdi.hdi_compound_mapping cm
  ON cm.source = 'np'
 AND cm.source_id = s.np_id
JOIN hdi.hdi_stg_final_target_chembl_data tmap
  ON tmap.target_id = s.target_id
JOIN hdi.hdi_target t
  ON t.source = 'chembl'
 AND t.source_target_id = tmap.chembl_id
GROUP BY
    cm.compound_key,
    t.target_key;

CREATE TABLE hdi.hdi_compound_metabolism (
    substrate_compound_key   VARCHAR(40) 
        REFERENCES hdi.hdi_compound (compound_key),
    metabolite_compound_key  VARCHAR(40)
        REFERENCES hdi.hdi_compound (compound_key),
    metabolism_target_key    BIGINT
        REFERENCES hdi.hdi_target (target_key),
    enzyme_name              TEXT,
    organism                 TEXT,
    tax_id                   INT,
    source                   TEXT NOT NULL DEFAULT 'chembl'
);
INSERT INTO hdi.hdi_compound_metabolism (
    substrate_compound_key,
    metabolite_compound_key,
    metabolism_target_key,
    enzyme_name,
    organism,
    tax_id,
    source
)
SELECT DISTINCT
    sub_map.compound_key       AS substrate_compound_key,
    met_map.compound_key       AS metabolite_compound_key,
    ht.target_key              AS metabolism_target_key,
    m.enzyme_name,
    m.organism,
    m.tax_id,
    'chembl'                   AS source

FROM public.metabolism m
-- substrate
JOIN public.compound_records sub_cr
  ON m.substrate_record_id = sub_cr.record_id
LEFT JOIN hdi.hdi_compound_mapping sub_map
  ON sub_cr.molregno = sub_map.source_sub_id::INT
 AND sub_map.source = 'chembl'
-- metabolite (LEFT เพราะอาจไม่มี)
LEFT JOIN public.compound_records met_cr
  ON m.metabolite_record_id = met_cr.record_id
LEFT JOIN hdi.hdi_compound_mapping met_map
  ON met_cr.molregno = met_map.source_sub_id::INT
 AND met_map.source = 'chembl'
-- metabolism target (LEFT เพราะ enzyme_tid อาจ NULL)
LEFT JOIN hdi.hdi_target ht
  ON m.enzyme_tid = ht.tid;

CREATE TABLE hdi.stg_drugbank_enzyme_action (
    drugbank_id   TEXT,
    drug_name     TEXT,
    chembl_id     TEXT,
    enzyme_name   TEXT,
    uniprot_id    TEXT,
    gene_name     TEXT,
    action        TEXT,
    organism      TEXT,
    source        TEXT
);
SELECT action, COUNT(*)
FROM hdi.stg_drugbank_enzyme_action
GROUP BY action
ORDER BY COUNT(*) DESC;

CREATE TABLE hdi.stg_drugbank_metabolism AS
SELECT *
FROM hdi.stg_drugbank_enzyme_action
WHERE action IN (
    'substrate',
    'metabolizer',
    'product of',
    'cleavage'
);
DELETE FROM hdi.stg_drugbank_metabolism
WHERE chembl_id IS NULL;
CREATE TABLE hdi.stg_drugbank_metabolism_mapped_drug AS
SELECT
    s.*,
    cm.compound_key
FROM hdi.stg_drugbank_metabolism s
JOIN hdi.hdi_compound_mapping cm
  ON cm.source = 'chembl'
 AND cm.source_id = s.chembl_id;
CREATE TABLE hdi.stg_uniprot_to_tid AS
SELECT DISTINCT
    cs.accession        AS uniprot_id,
    tc.tid              AS tid,
    cs.organism,
    cs.tax_id
FROM public.component_sequences cs
JOIN public.target_components tc
  ON cs.component_id = tc.component_id
WHERE cs.component_type = 'PROTEIN';

CREATE TABLE hdi.stg_uniprot_single_protein_tid AS
SELECT
    u.uniprot_id,
    u.tid,
    u.organism,
    u.tax_id
FROM hdi.stg_uniprot_to_tid u
JOIN (
    SELECT
        tc.tid
    FROM public.target_components tc
    JOIN public.component_sequences cs
      ON tc.component_id = cs.component_id
    WHERE cs.component_type = 'PROTEIN'
    GROUP BY tc.tid
    HAVING COUNT(DISTINCT tc.component_id) = 1
) single_protein
  ON u.tid = single_protein.tid;

SELECT
    uniprot_id,
    COUNT(DISTINCT tid)
FROM hdi.stg_uniprot_to_tid
GROUP BY uniprot_id
HAVING COUNT(DISTINCT tid) > 1;

CREATE TABLE hdi.stg_uniprot_to_target_key AS
SELECT
    s.uniprot_id,
    ht.target_key,
    ht.source_target_id AS chembl_target_id,
    ht.pref_name        AS target_name,
    ht.organism,
    ht.tax_id
FROM hdi.stg_uniprot_single_protein_tid s
JOIN hdi.hdi_target ht
  ON ht.tid = s.tid
 AND ht.source = 'chembl'
WHERE ht.is_protein = TRUE;
INSERT INTO hdi.hdi_compound_metabolism (
    substrate_compound_key,
    metabolite_compound_key,
    metabolism_target_key,
    enzyme_name,
    organism,
    tax_id,
    source
)
SELECT DISTINCT
    d.compound_key      AS substrate_compound_key,
    NULL                AS metabolite_compound_key,
    u.target_key       AS metabolism_target_key,
    d.enzyme_name,
    u.organism,
    u.tax_id,
    'drugbank'         AS source
FROM hdi.stg_drugbank_metabolism_mapped_drug d
JOIN hdi.stg_uniprot_to_target_key u
  ON d.uniprot_id = u.uniprot_id
WHERE d.action IN ('substrate', 'metabolizer');
INSERT INTO hdi.hdi_compound_metabolism (
    substrate_compound_key,
    metabolite_compound_key,
    metabolism_target_key,
    enzyme_name,
    organism,
    tax_id,
    source
)
SELECT DISTINCT
    NULL                AS substrate_compound_key,
    d.compound_key      AS metabolite_compound_key,
    u.target_key       AS metabolism_target_key,
    d.enzyme_name,
    u.organism,
    u.tax_id,
    'drugbank'         AS source
FROM hdi.stg_drugbank_metabolism_mapped_drug d
JOIN hdi.stg_uniprot_to_target_key u
  ON d.uniprot_id = u.uniprot_id
WHERE d.action IN ('product of', 'cleavage'); 

CREATE TABLE hdi.stg_drugbank_enzyme_effect AS
SELECT *
FROM hdi.stg_drugbank_enzyme_action
WHERE action IN (
    'inhibitor',
    'inducer',
    'activator',
    'downregulator',
	'inactivator'
)
AND chembl_id IS NOT NULL;
CREATE TABLE hdi.stg_drugbank_enzyme_effect_mapped_drug AS
SELECT
    s.*,
    cm.compound_key
FROM hdi.stg_drugbank_enzyme_effect s
JOIN hdi.hdi_compound_mapping cm
  ON cm.source = 'chembl'
 AND cm.source_id = s.chembl_id;

CREATE TABLE hdi.hdi_compound_enzyme_effect (
    compound_key      VARCHAR(40)
        REFERENCES hdi.hdi_compound (compound_key),
    enzyme_target_key BIGINT
        REFERENCES hdi.hdi_target (target_key),
    effect_type       TEXT
        CHECK (effect_type IN (
            'inhibitor',
            'inducer',
            'activator',
            'downregulator',
			'inactivator'
        )),
    enzyme_name       TEXT,
    gene_name         TEXT,
    organism          TEXT,
    tax_id            INT,
    source            TEXT NOT NULL DEFAULT 'drugbank',
    PRIMARY KEY (compound_key, enzyme_target_key, effect_type, source)
);

INSERT INTO hdi.hdi_compound_enzyme_effect (
    compound_key,
    enzyme_target_key,
    effect_type,
    enzyme_name,
    gene_name,
    organism,
    tax_id,
    source
)
SELECT DISTINCT
    d.compound_key,
    u.target_key            AS enzyme_target_key,

    CASE
        WHEN d.action IN ('inhibitor')      THEN 'inhibitor'
        WHEN d.action IN ('inducer')        THEN 'inducer'
        WHEN d.action IN ('activator')      THEN 'activator'
        WHEN d.action IN ('downregulator')  THEN 'downregulator'
		WHEN d.action IN ('inactivator')    THEN 'inactivator'
		
    END                     AS effect_type,

    d.enzyme_name,
    d.gene_name,
    d.organism,
    u.tax_id,
    'drugbank'              AS source

FROM hdi.stg_drugbank_enzyme_effect_mapped_drug d
JOIN hdi.stg_uniprot_to_target_key u
  ON u.uniprot_id = d.uniprot_id
WHERE d.action IN (
    'inhibitor',
    'inducer',
    'activator',
    'downregulator',
    'inactivator'
)
ON CONFLICT (compound_key, enzyme_target_key, effect_type, source)
DO NOTHING;
CREATE TABLE hdi.hdi_atc (
    atc_code        TEXT PRIMARY KEY,
    atc_level       INT,
    atc_name        TEXT
);
INSERT INTO hdi.hdi_atc (atc_code, atc_level, atc_name)
SELECT DISTINCT level1, 1, level1_description
FROM public.atc_classification
WHERE level1 IS NOT NULL

UNION
SELECT DISTINCT level2, 2, level2_description
FROM public.atc_classification
WHERE level2 IS NOT NULL

UNION
SELECT DISTINCT level3, 3, level3_description
FROM public.atc_classification
WHERE level3 IS NOT NULL

UNION
SELECT DISTINCT level4, 4, level4_description
FROM public.atc_classification
WHERE level4 IS NOT NULL

UNION
SELECT DISTINCT level5, 5, who_name
FROM public.atc_classification
WHERE level5 IS NOT NULL;

CREATE TABLE hdi.hdi_atc_hierarchy (
    child_atc   TEXT REFERENCES hdi.hdi_atc(atc_code),
    parent_atc  TEXT REFERENCES hdi.hdi_atc(atc_code),
    PRIMARY KEY (child_atc, parent_atc)
);
INSERT INTO hdi.hdi_atc_hierarchy (child_atc, parent_atc)
SELECT DISTINCT level5, level4
FROM public.atc_classification
WHERE level5 IS NOT NULL AND level4 IS NOT NULL

UNION
SELECT DISTINCT level4, level3
FROM public.atc_classification
WHERE level4 IS NOT NULL AND level3 IS NOT NULL

UNION
SELECT DISTINCT level3, level2
FROM public.atc_classification
WHERE level3 IS NOT NULL AND level2 IS NOT NULL

UNION
SELECT DISTINCT level2, level1
FROM public.atc_classification
WHERE level2 IS NOT NULL AND level1 IS NOT NULL;

CREATE TABLE hdi.hdi_compound_atc (
    compound_key VARCHAR(40)
        REFERENCES hdi.hdi_compound(compound_key),
    atc5_code    TEXT
        REFERENCES hdi.hdi_atc(atc_code),
    source       TEXT DEFAULT 'chembl',
    PRIMARY KEY (compound_key, atc5_code, source)
);
INSERT INTO hdi.hdi_compound_atc (
    compound_key,
    atc5_code,
    source
)
SELECT DISTINCT
    cm.compound_key,
    m.level5,
    'chembl'
FROM public.molecule_atc_classification m
JOIN hdi.hdi_compound_mapping cm
  ON cm.source = 'chembl'
 AND cm.source_sub_id = m.molregno::TEXT;

CREATE TABLE hdi.stg_species_info_selected (
    org_id TEXT,
    org_name TEXT,
    org_tax_id TEXT,
    species_tax_id TEXT,
    species_name TEXT,
    genus_tax_id TEXT,
    genus_name TEXT,
    family_tax_id TEXT,
    family_name TEXT,
    kingdom_name TEXT,
    kingdom_tax_id TEXT,
    superkingdom_name TEXT,
    superkingdom_tax_id TEXT
);

CREATE TABLE hdi.stg_species_pair_selected (
    org_id TEXT,
    np_id TEXT
);

CREATE TABLE hdi.P_organism ( 
	org_id TEXT PRIMARY KEY, 
	org_name TEXT, 
	org_tax_id TEXT );

INSERT INTO hdi.P_organism (org_id, org_name, org_tax_id)
SELECT DISTINCT
    org_id,
    org_name,
    org_tax_id
FROM hdi.stg_species_info_selected
WHERE org_id IS NOT NULL
  AND org_id <> ''
  AND org_id <> 'n.a.';
  
CREATE TABLE hdi.P_species (
    species_tax_id TEXT PRIMARY KEY,
    species_name   TEXT
);
INSERT INTO hdi.P_species (species_tax_id, species_name)
SELECT
    species_tax_id,
    MIN(species_name) AS species_name
FROM hdi.stg_species_info_selected
WHERE species_tax_id IS NOT NULL
  AND species_tax_id <> ''
  AND species_tax_id <> 'n.a.'
GROUP BY species_tax_id;

CREATE TABLE hdi.P_genus (
    genus_tax_id TEXT PRIMARY KEY,
    genus_name   TEXT
);
INSERT INTO hdi.P_genus (genus_tax_id, genus_name)
SELECT
    genus_tax_id,
    MIN(genus_name) AS genus_name
FROM hdi.stg_species_info_selected
WHERE genus_tax_id IS NOT NULL
  AND genus_tax_id <> ''
  AND genus_tax_id <> 'n.a.'
GROUP BY genus_tax_id;

CREATE TABLE hdi.P_family (
    family_tax_id TEXT PRIMARY KEY,
    family_name   TEXT
);
INSERT INTO P_family (family_tax_id, family_name)
SELECT
    family_tax_id,
    MIN(family_name) AS family_name
From hdi.stg_species_info_selected
WHERE family_tax_id IS NOT NULL
  AND family_tax_id <> ''
  AND family_tax_id <> 'n.a.'
GROUP BY family_tax_id;

CREATE TABLE P_kingdom (
    kingdom_tax_id TEXT PRIMARY KEY,
    kingdom_name   TEXT
);
INSERT INTO P_kingdom (kingdom_tax_id, kingdom_name)
SELECT
    kingdom_tax_id,
    MIN(kingdom_name) AS kingdom_name
FROM hdi.stg_species_info_selected
WHERE kingdom_tax_id IS NOT NULL
  AND kingdom_tax_id <> ''
  AND kingdom_tax_id <> 'n.a.'
GROUP BY kingdom_tax_id;

CREATE TABLE P_superkingdom (
    superkingdom_tax_id TEXT PRIMARY KEY,
    superkingdom_name   TEXT
);
INSERT INTO P_superkingdom (superkingdom_tax_id, superkingdom_name)
SELECT
    superkingdom_tax_id,
    MIN(superkingdom_name) AS superkingdom_name
FROM hdi.stg_species_info_selected
WHERE superkingdom_tax_id IS NOT NULL
  AND superkingdom_tax_id <> ''
  AND superkingdom_tax_id <> 'n.a.'
GROUP BY superkingdom_tax_id;

CREATE TABLE hdi.R_organism_species (
    org_id         TEXT,
    species_tax_id TEXT,
    PRIMARY KEY (org_id, species_tax_id),
    FOREIGN KEY (org_id) REFERENCES hdi.P_organism(org_id),
    FOREIGN KEY (species_tax_id) REFERENCES hdi.P_species(species_tax_id)
);
INSERT INTO hdi.R_organism_species (org_id, species_tax_id)
SELECT DISTINCT
    org_id,
    species_tax_id
FROM hdi.stg_species_info_selected
WHERE org_id IS NOT NULL
  AND species_tax_id IS NOT NULL
  AND org_id <> ''
  AND species_tax_id <> ''
  AND org_id <> 'n.a.'
  AND species_tax_id <> 'n.a.';
CREATE TABLE hdi.R_species_genus (
    species_tax_id TEXT,
    genus_tax_id   TEXT,
    PRIMARY KEY (species_tax_id, genus_tax_id),
    FOREIGN KEY (species_tax_id) REFERENCES hdi.P_species(species_tax_id),
    FOREIGN KEY (genus_tax_id) REFERENCES hdi.P_genus(genus_tax_id)
);
INSERT INTO hdi.R_species_genus (species_tax_id, genus_tax_id)
SELECT DISTINCT
    species_tax_id,
    genus_tax_id
FROM hdi.stg_species_info_selected
WHERE species_tax_id IS NOT NULL
  AND genus_tax_id IS NOT NULL
  AND species_tax_id <> 'n.a.'
  AND genus_tax_id <> 'n.a.';
CREATE TABLE hdi.R_genus_family (
    genus_tax_id  TEXT,
    family_tax_id TEXT,
    PRIMARY KEY (genus_tax_id, family_tax_id),
    FOREIGN KEY (genus_tax_id) REFERENCES hdi.P_genus(genus_tax_id),
    FOREIGN KEY (family_tax_id) REFERENCES hdi.P_family(family_tax_id)
);
INSERT INTO hdi.R_genus_family (genus_tax_id, family_tax_id)
SELECT DISTINCT
    genus_tax_id,
    family_tax_id
FROM hdi.stg_species_info_selected
WHERE genus_tax_id IS NOT NULL
  AND family_tax_id IS NOT NULL
  AND genus_tax_id <> 'n.a.'
  AND family_tax_id <> 'n.a.';
 
CREATE TABLE hdi.R_family_kingdom (
    family_tax_id  TEXT,
    kingdom_tax_id TEXT,
    PRIMARY KEY (family_tax_id, kingdom_tax_id),
    FOREIGN KEY (family_tax_id) REFERENCES hdi.P_family(family_tax_id),
    FOREIGN KEY (kingdom_tax_id) REFERENCES hdi.P_kingdom(kingdom_tax_id)
);
INSERT INTO hdi.R_family_kingdom (family_tax_id, kingdom_tax_id)
SELECT DISTINCT
    family_tax_id,
    kingdom_tax_id
FROM hdi.stg_species_info_selected
WHERE family_tax_id IS NOT NULL
  AND kingdom_tax_id IS NOT NULL
  AND family_tax_id <> 'n.a.'
  AND kingdom_tax_id <> 'n.a.';
CREATE TABLE hdi.R_kingdom_superkingdom (
    kingdom_tax_id      TEXT,
    superkingdom_tax_id TEXT,
    PRIMARY KEY (kingdom_tax_id, superkingdom_tax_id),
    FOREIGN KEY (kingdom_tax_id) REFERENCES hdi.P_kingdom(kingdom_tax_id),
    FOREIGN KEY (superkingdom_tax_id) REFERENCES hdi.P_superkingdom(superkingdom_tax_id)
);
INSERT INTO hdi.R_kingdom_superkingdom (kingdom_tax_id, superkingdom_tax_id)
SELECT DISTINCT
    kingdom_tax_id,
    superkingdom_tax_id
FROM hdi.stg_species_info_selected
WHERE kingdom_tax_id IS NOT NULL
  AND superkingdom_tax_id IS NOT NULL
  AND kingdom_tax_id <> 'n.a.'
  AND superkingdom_tax_id <> 'n.a.';

CREATE TABLE hdi.R_organism_compound (
    org_id       TEXT
        REFERENCES hdi.P_organism (org_id),
    compound_key TEXT
        REFERENCES hdi.hdi_compound (compound_key),
    PRIMARY KEY (org_id, compound_key)
); 
INSERT INTO hdi.R_organism_compound (
    org_id,
    compound_key
)
SELECT DISTINCT
    s.org_id,
    m.compound_key
FROM hdi.stg_species_pair_selected s
JOIN hdi.hdi_compound_mapping m
  ON m.source = 'np'
 AND m.source_id = s.np_id
JOIN hdi.P_organism o
  ON o.org_id = s.org_id
ON CONFLICT DO NOTHING;

