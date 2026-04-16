CREATE TABLE hdi.ddi_chembl_pairs (
    chembl_id_1 TEXT NOT NULL,
    chembl_id_2 TEXT NOT NULL,
    -- ห้าม self interaction
    CONSTRAINT no_self_interaction 
        CHECK (chembl_id_1 <> chembl_id_2),
    -- บังคับเรียงลำดับ
    CONSTRAINT enforce_order
        CHECK (chembl_id_1 < chembl_id_2),
    -- กันซ้ำ
    PRIMARY KEY (chembl_id_1, chembl_id_2)
);
 
CREATE TABLE hdi.ddi_compound_pairs (
    compound_key_1 VARCHAR(40) NOT NULL,
    compound_key_2 VARCHAR(40) NOT NULL,
    CONSTRAINT no_self_loop 
        CHECK (compound_key_1 <> compound_key_2),
    CONSTRAINT enforce_order
        CHECK (compound_key_1 < compound_key_2),
    PRIMARY KEY (compound_key_1, compound_key_2)
);
INSERT INTO hdi.ddi_compound_pairs
SELECT DISTINCT
    LEAST(m1.compound_key, m2.compound_key),
    GREATEST(m1.compound_key, m2.compound_key)
FROM hdi.ddi_chembl_pairs d
JOIN hdi.hdi_compound_mapping m1
    ON d.chembl_id_1 = m1.source_id
    AND m1.source = 'chembl'
JOIN hdi.hdi_compound_mapping m2
    ON d.chembl_id_2 = m2.source_id
    AND m2.source = 'chembl';

