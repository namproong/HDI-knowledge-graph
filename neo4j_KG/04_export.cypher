CALL apoc.export.csv.query(
"
MATCH (h)-[r]->(t)


WITH h, r, t


WITH
CASE
    WHEN h:Compound THEN toString(h.id)
    WHEN h:Target THEN toString(h.id)
    WHEN h:ATC1 THEN toString(h.id)
    WHEN h:ATC2 THEN toString(h.id)
    WHEN h:ATC3 THEN toString(h.id)
    WHEN h:ATC4 THEN toString(h.id)
    WHEN h:ATC5 THEN toString(h.id)
    WHEN h:ProteinClass THEN 'PC_' + toString(h.protein_class_id)
    WHEN h:Organism THEN toString(h.org_id)
    WHEN h:Species THEN 'SPC_' + toString(h.species_tax_id)
    WHEN h:Genus THEN 'GEN_' + toString(h.genus_tax_id)
    WHEN h:Family THEN 'FAM_' + toString(h.family_tax_id)
    WHEN h:Kingdom THEN 'KING_' + toString(h.kingdom_tax_id)
    WHEN h:Superkingdom THEN 'SKING_' + toString(h.superkingdom_tax_id)
END AS head_id,


type(r) AS relation,


CASE
    WHEN t:Compound THEN toString(t.id)
    WHEN t:Target THEN toString(t.id)
    WHEN t:ATC1 THEN toString(t.id)
    WHEN t:ATC2 THEN toString(t.id)
    WHEN t:ATC3 THEN toString(t.id)
    WHEN t:ATC4 THEN toString(t.id)
    WHEN t:ATC5 THEN toString(t.id)
    WHEN t:ProteinClass THEN 'PC_' + toString(t.protein_class_id)
    WHEN t:Organism THEN toString(t.org_id)
    WHEN t:Species THEN 'SPC_' + toString(t.species_tax_id)
    WHEN t:Genus THEN 'GEN_' + toString(t.genus_tax_id)
    WHEN t:Family THEN 'FAM_' + toString(t.family_tax_id)
    WHEN t:Kingdom THEN 'KING_' + toString(t.kingdom_tax_id)
    WHEN t:Superkingdom THEN 'SKING_' + toString(t.superkingdom_tax_id)
END AS tail_id


RETURN head_id, relation, tail_id
",
"HDI_triples_finalV2forP.csv",
