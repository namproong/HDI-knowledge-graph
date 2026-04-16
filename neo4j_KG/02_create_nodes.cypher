LOAD CSV WITH HEADERS FROM 'file:///hdi_compound.csv' AS row
CALL {
  WITH row
  MERGE (c:Compound {id: row.compound_key})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///hdi_target.csv' AS row
CALL (row) {
  MERGE (t:Target {id: row.target_key})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc.csv' AS row
CALL (row) {
  WITH row
  WHERE row.atc_level = '5'
  MERGE (:ATC5 {id: row.atc_code})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc.csv' AS row
CALL (row) {
  WITH row
  WHERE row.atc_level = '4'
  MERGE (:ATC4 {id: row.atc_code})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc.csv' AS row
CALL (row) {
  WITH row
  WHERE row.atc_level = '3'
  MERGE (:ATC3 {id: row.atc_code})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc.csv' AS row
CALL (row) {
  WITH row
  WHERE row.atc_level = '2'
  MERGE (:ATC2{id: row.atc_code})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc.csv' AS row
CALL (row) {
  WITH row
  WHERE row.atc_level = '1'
  MERGE (:ATC1 {id: row.atc_code})
} IN TRANSACTIONS OF 1000 ROWS;

LOAD CSV WITH HEADERS
FROM 'file:///hdi_protein_class_hierarchy.csv'
AS row
MERGE (:ProteinClass {protein_class_id: row.protein_class_id})
MERGE (:ProteinClass {protein_class_id: row.parent_id});

LOAD CSV WITH HEADERS
FROM 'file:///p_superkingdom.csv'
AS row
MERGE (:Superkingdom {superkingdom_tax_id: row.superkingdom_tax_id, name: row.superkingdom_name});

LOAD CSV WITH HEADERS
FROM 'file:///p_kingdom.csv'
AS row
MERGE (:Kingdom {kingdom_tax_id: row.kingdom_tax_id, name: row.kingdom_name});

LOAD CSV WITH HEADERS FROM 'file:///p_family.csv' AS row
MERGE (:Family {family_tax_id: row.family_tax_id, name: row.family_name});

LOAD CSV WITH HEADERS FROM 'file:///p_genus.csv' AS row
MERGE (:Genus {genus_tax_id: row.genus_tax_id, name: row.genus_name});

LOAD CSV WITH HEADERS FROM 'file:///p_species.csv' AS row
MERGE (:Species {species_tax_id: row.species_tax_id, name: row.species_name});

LOAD CSV WITH HEADERS FROM 'file:///p_organism.csv' AS row
MERGE (:Organism {org_id: row.org_id, name: row.org_name});
