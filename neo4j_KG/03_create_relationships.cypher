LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'NEGATIVE_MODULATOR'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:NEGATIVELY_MODULATES]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'POSITIVE_MODULATOR'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:POSITIVELY_MODULATES]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'UNSPECIFIED_MODULATION'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:MODULATES]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'BINDING_ONLY'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:BINDS]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'ENZYMATIC_ACTIVITY'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:AFFECTS_ENZYMATIC_FUNCTION]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'STRUCTURAL_MODULATION'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:STRUCTURALLY_MODULATES]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_compound_target_pd.csv' AS row
WITH row
WHERE row.parent_action_type = 'SUBSTRATE_TRANSPORT'
MATCH (c:Compound {id: row.compound_key})
MATCH (t:Target {id: row.target_key})
MERGE (c)-[:IS_TRANSPORT_SUBSTRATE_OF]->(t);

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc_hierarchy.csv' AS row
WITH row
WHERE size(row.child_atc) = 7 AND size(row.parent_atc) = 5
MATCH (c:ATC5 {id: row.child_atc})
MATCH (p:ATC4 {id: row.parent_atc})
MERGE (c)-[:SUBCLASS_OF_ATC4]->(p);

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc_hierarchy.csv' AS row
WITH row
WHERE size(row.child_atc) = 5 AND size(row.parent_atc) = 4
MATCH (c:ATC4 {id: row.child_atc})
MATCH (p:ATC3 {id: row.parent_atc})
MERGE (c)-[:SUBCLASS_OF_ATC3]->(p);

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc_hierarchy.csv' AS row
WITH row
WHERE size(row.child_atc) = 4 AND size(row.parent_atc) = 3
MATCH (c:ATC3 {id: row.child_atc})
MATCH (p:ATC2 {id: row.parent_atc})
MERGE (c)-[:SUBCLASS_OF_ATC2]->(p);

LOAD CSV WITH HEADERS FROM 'file:///hdi_atc_hierarchy.csv' AS row
WITH row
WHERE size(row.child_atc) = 3 AND size(row.parent_atc) = 1
MATCH (c:ATC2 {id: row.child_atc})
MATCH (p:ATC1 {id: row.parent_atc})
MERGE (c)-[:SUBCLASS_OF_ATC1]->(p);

LOAD CSV WITH HEADERS
FROM 'file:///hdi_compound_atc.csv'
AS row
MATCH (c:Compound {id: row.compound_key})
MATCH (a:ATC5 {id: row.atc5_code})
MERGE (c)-[:HAS_ATC5]->(a);

LOAD CSV WITH HEADERS
FROM 'file:///hdi_target_protein_class.csv'
AS row
MATCH (t:Target {id: row.target_key})
MATCH (pc:ProteinClass {protein_class_id: row.protein_class_id})
MERGE (t)-[:HAS_PROTEIN_CLASS]->(pc);

LOAD CSV WITH HEADERS
FROM 'file:///hdi_protein_class_hierarchy.csv'
AS row
MATCH (child:ProteinClass {protein_class_id: row.protein_class_id})
MATCH (parent:ProteinClass {protein_class_id: row.parent_id})
MERGE (child)-[:IS_SUBCLASS_OF_PC]->(parent);


LOAD CSV WITH HEADERS
FROM 'file:///r_species_genus.csv'
AS row
MATCH (s:Species {species_tax_id: row.species_tax_id})
MATCH (g:Genus {genus_tax_id: row.genus_tax_id})
MERGE (s)-[:BELONGS_TO_GENUS]->(g); 

LOAD CSV WITH HEADERS
FROM 'file:///r_genus_family.csv'
AS row
MATCH (g:Genus {genus_tax_id: row.genus_tax_id})
MATCH (f:Family {family_tax_id: row.family_tax_id})
MERGE (g)-[:BELONGS_TO_FAMILY]->(f);

LOAD CSV WITH HEADERS
FROM 'file:///r_family_kingdom.csv'
AS row
MATCH (f:Family {family_tax_id: row.family_tax_id})
MATCH (k:Kingdom {kingdom_tax_id: row.kingdom_tax_id})
MERGE (f)-[:BELONGS_TO_KINGDOM]->(k);

LOAD CSV WITH HEADERS
FROM 'file:///r_kingdom_superkingdom.csv'
AS row
MATCH (k:Kingdom {kingdom_tax_id: row.kingdom_tax_id})
MATCH (sk:Superkingdom {superkingdom_tax_id: row.superkingdom_tax_id})
MERGE (k)-[:BELONGS_TO_SUPERKINGDOM]->(sk);

CALL apoc.periodic.iterate(
  "LOAD CSV WITH HEADERS
   FROM 'file:///r_organism_compound.csv'
   AS row RETURN row",
  "MATCH (c:Compound {id: row.compound_key})
   MATCH (o:Organism {org_id: row.org_id})
   MERGE (c)-[:FOUND_IN_ORGANISM]->(o)",
  {batchSize: 2000, parallel: false}
);

LOAD CSV WITH HEADERS
FROM 'file:///r_organism_species.csv'
AS row
MATCH (o:Organism {org_id: row.org_id})
MATCH (s:Species {species_tax_id: row.species_tax_id})
CREATE (o)-[:HAS_SPECIES]->(s);

CALL apoc.periodic.iterate(
  "
  LOAD CSV WITH HEADERS
  FROM 'file:///hdi_compound_enzyme_effect.csv'
  AS row
  RETURN row
  ",
  "
  MATCH (c:Compound {id: row.compound_key})
  MATCH (t:Target   {id: row.enzyme_target_key})

  WITH c, t,
       trim(toUpper(row.effect_type)) AS effect

  FOREACH (_ IN CASE WHEN effect = 'INHIBITOR' THEN [1] ELSE [] END |
    MERGE (c)-[:INHIBITS]->(t)
  )
  FOREACH (_ IN CASE WHEN effect = 'INDUCER' THEN [1] ELSE [] END |
    MERGE (c)-[:INDUCES]->(t)
  )
  FOREACH (_ IN CASE WHEN effect = 'ACTIVATOR' THEN [1] ELSE [] END |
    MERGE (c)-[:ACTIVATES]->(t)
  )
  FOREACH (_ IN CASE WHEN effect = 'DOWNREGULATOR' THEN [1] ELSE [] END |
    MERGE (c)-[:DOWNREGULATES]->(t)
  )
  ",
  {batchSize: 2000, parallel: false}
);

CALL apoc.periodic.iterate(
  "
  LOAD CSV WITH HEADERS
  FROM 'file:///hdi_compound_metabolism.csv'
  AS row
  RETURN row
  ",
  "
  MATCH (s:Compound {id: row.substrate_compound_key})
  MATCH (m:Compound {id: row.metabolite_compound_key})
  MERGE (s)-[:METABOLIZED_TO]->(m)
  ",
  {batchSize: 2000, parallel: false}
);

CALL apoc.periodic.iterate(
  "
  LOAD CSV WITH HEADERS
  FROM 'file:///hdi_compound_metabolism.csv'
  AS row
  WITH row
  WHERE row.metabolism_target_key IS NOT NULL
  RETURN row
  ",
  "
  MATCH (s:Compound {id: row.substrate_compound_key})
  MATCH (t:Target   {id: row.metabolism_target_key})
  MERGE (s)-[:METABOLIZED_BY]->(t)
  ",
  {batchSize: 2000, parallel: false}
);

CALL apoc.periodic.iterate(
  "
  LOAD CSV WITH HEADERS
  FROM 'file:///hdi_compound_target_assoc.csv'
  AS row
  RETURN row
  ",
  "
  MATCH (c:Compound {id: row.compound_key})
  MATCH (t:Target{id: row.target_key})

  MERGE (c)-[r:INTERACTS_WITH]->(t)
  ON CREATE SET
    r.source = row.source,
    r.evidence_count = toInteger(row.evidence_count)
  ON MATCH SET
    r.evidence_count =
      coalesce(r.evidence_count,0) + toInteger(row.evidence_count)
  ",
  {batchSize: 2000, parallel: false}
);