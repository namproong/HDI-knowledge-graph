CREATE CONSTRAINT compound_id_unique
FOR (c:Compound)
REQUIRE c.id IS UNIQUE;

CREATE CONSTRAINT target_id_unique
FOR (t:Target)
REQUIRE t.id IS UNIQUE;

CREATE INDEX organism_id_index IF NOT EXISTS
FOR (o:Organism)
ON (o.id);