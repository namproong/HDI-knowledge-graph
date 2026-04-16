SELECT kingdom_name,
       COUNT(*) AS count_rows
FROM hdi.stg_species_info_selected
GROUP BY kingdom_name
ORDER BY count_rows DESC;

SELECT COUNT(DISTINCT roc.compound_key) AS plant_compound_count
FROM hdi.r_organism_compound roc
JOIN hdi.stg_species_info_selected s
  ON roc.org_id = s.org_id
WHERE s.kingdom_name = 'Viridiplantae'; 116028 

SELECT COUNT(DISTINCT compound_key) AS total_with_organism
FROM hdi.r_organism_compound; 201827

SELECT s.kingdom_name,
       COUNT(DISTINCT roc.compound_key) AS compound_count
FROM hdi.r_organism_compound roc
JOIN hdi.stg_species_info_selected s
  ON roc.org_id = s.org_id
GROUP BY s.kingdom_name
ORDER BY compound_count DESC;
"Viridiplantae"	116028
"Metazoa"	47381
"Fungi"	47212
"n.a."	44060
"Bacillati"	11497
"Pseudomonadati"	1865
"Methanobacteriati"	48
"Thermotogati"	4
"Thermoproteati"	2
 

SELECT roc.compound_key,
       COUNT(DISTINCT s.kingdom_name) AS kingdom_count
FROM hdi.r_organism_compound roc
JOIN hdi.stg_species_info_selected s
  ON roc.org_id = s.org_id
GROUP BY roc.compound_key
HAVING COUNT(DISTINCT s.kingdom_name) > 1; มี 34666 results 

SELECT COUNT(DISTINCT roc.compound_key)
FROM hdi.r_organism_compound roc
JOIN hdi.stg_species_info_selected s
  ON roc.org_id = s.org_id
JOIN hdi.hdi_compound m
  ON roc.compound_key = m.compound_key
WHERE s.kingdom_name = 'Viridiplantae'; 116028

CREATE TABLE hdi.plant_compound AS
SELECT DISTINCT roc.compound_key
FROM hdi.r_organism_compound roc
JOIN hdi.stg_species_info_selected s
  ON roc.org_id = s.org_id
WHERE s.kingdom_name = 'Viridiplantae';

SELECT COUNT(*) FROM hdi.plant_compound; 116028

SELECT COUNT(DISTINCT compound_key) AS plant_in_ddi
FROM (
    SELECT d.compound_key_1 AS compound_key
    FROM hdi.ddi_compound_pairs d
    JOIN hdi.plant_compound p
      ON d.compound_key_1 = p.compound_key

    UNION

    SELECT d.compound_key_2
    FROM hdi.ddi_compound_pairs d
    JOIN hdi.plant_compound p
      ON d.compound_key_2 = p.compound_key
) t; 528

SELECT
    CASE
        WHEN p1.compound_key IS NOT NULL
         AND p2.compound_key IS NOT NULL THEN 'plant-plant'
        WHEN p1.compound_key IS NOT NULL
         OR p2.compound_key IS NOT NULL THEN 'plant-nonplant'
        ELSE 'nonplant-nonplant'
    END AS pair_type,
    COUNT(*) AS pair_count
FROM hdi.ddi_compound_pairs d
LEFT JOIN hdi.plant_compound p1
  ON d.compound_key_1 = p1.compound_key
LEFT JOIN hdi.plant_compound p2
  ON d.compound_key_2 = p2.compound_key
GROUP BY pair_type;
"pair_type"	"pair_count"
"nonplant-nonplant"	761048
"plant-plant"	25150
"plant-nonplant"	273218

CREATE TABLE hdi.species_compound_notF AS
SELECT DISTINCT
    ps.species_tax_id,
    ps.species_name,
    roc.compound_key
FROM hdi.r_organism_compound roc
JOIN hdi.r_organism_species ros
    ON roc.org_id = ros.org_id
JOIN hdi.p_species ps
    ON ros.species_tax_id = ps.species_tax_id
ORDER BY ps.species_tax_id, roc.compound_key;