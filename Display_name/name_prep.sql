-- 1. สร้างตารางเก็บข้อมูล Natural Products เบื้องต้น
CREATE TABLE hdi.natural_products_selected ( 
    np_id TEXT, 
    pref_name TEXT, 
    iupac_name TEXT 
);

-- 2. สร้างตารางชื่อสารประกอบ โดยคลีนข้อมูลชื่อที่ใช้ไม่ได้ (เช่น n.a., no data) ให้เป็น NULL
CREATE TABLE hdi.npass_compound_names AS 
SELECT 
    m.compound_key, 
    n.np_id, 
    CASE 
        WHEN LOWER(TRIM(n.pref_name)) IN ( 'n.a.', 'na', 'no data', 'not named' ) THEN NULL 
        ELSE NULLIF(TRIM(n.pref_name), '') 
    END AS pref_name, 
    CASE 
        WHEN LOWER(TRIM(n.iupac_name)) IN ( 'n.a.', 'na', 'no data', 'not named' ) THEN NULL 
        ELSE NULLIF(TRIM(n.iupac_name), '') 
    END AS iupac_name 
FROM hdi.natural_products_selected n 
JOIN hdi.hdi_compound_mapping m 
    ON m.source = 'np' AND m.source_id = n.np_id; 

-- 3. เช็คจำนวน pref_name ที่ซ้ำกัน
SELECT pref_name, COUNT(*) 
FROM hdi.npass_compound_names 
WHERE pref_name IS NOT NULL 
GROUP BY pref_name 
HAVING COUNT(*) > 1 
ORDER BY COUNT(*) DESC; 
-- ผลลัพธ์: 1978 แถว 

-- 4. เช็คจำนวน iupac_name ที่ซ้ำกัน
SELECT iupac_name, COUNT(*) 
FROM hdi.npass_compound_names 
WHERE iupac_name IS NOT NULL 
GROUP BY iupac_name 
HAVING COUNT(*) > 1 
ORDER BY COUNT(*) DESC; 
-- ผลลัพธ์: 54 

-- 5. เช็คจำนวน pref_name ที่ซ้ำกัน โดยนับแยกตาม compound_key
SELECT pref_name, COUNT(DISTINCT compound_key) 
FROM hdi.npass_compound_names 
WHERE pref_name IS NOT NULL 
GROUP BY pref_name 
HAVING COUNT(DISTINCT compound_key) > 1 
ORDER BY COUNT(*) DESC; 
-- ผลลัพธ์: 1967 แถว 

-- 6. เช็คจำนวนรายการที่ไม่มีทั้ง pref_name และ iupac_name
SELECT COUNT(*) 
FROM hdi.npass_compound_names 
WHERE pref_name IS NULL AND iupac_name IS NULL; 
-- ผลลัพธ์: 3 

-- 7. ดูข้อมูลรายการที่ไม่มีทั้ง pref_name และ iupac_name
SELECT * FROM hdi.npass_compound_names 
WHERE pref_name IS NULL AND iupac_name IS NULL; 

-- 8. สร้างตารางสำหรับชื่อที่จะแสดงผล โดยเรียงลำดับความสำคัญ (pref_name -> iupac_name -> compound_key)
CREATE TABLE hdi.npass_display_name AS 
SELECT DISTINCT ON (compound_key) 
    compound_key, 
    np_id, 
    pref_name, 
    iupac_name, 
    COALESCE(pref_name, iupac_name, compound_key) AS display_name, 
    CASE 
        WHEN pref_name IS NOT NULL THEN 'pref_name' 
        WHEN iupac_name IS NOT NULL THEN 'iupac_name' 
        ELSE 'compound_key_fallback' 
    END AS name_source 
FROM hdi.npass_compound_names 
ORDER BY compound_key, pref_name IS NULL; 
-- ผลลัพธ์: 201862 

-- 9. เช็คว่ามีรายการไหนที่ไม่มี display_name หรือไม่
SELECT COUNT(*) 
FROM hdi.npass_display_name 
WHERE display_name IS NULL; 
-- ผลลัพธ์: 0 

-- 10. เช็คจำนวน compound_key ที่มาจาก chembl
SELECT COUNT(DISTINCT compound_key) 
FROM hdi.hdi_compound_name 
WHERE source = 'chembl'; 
-- ผลลัพธ์: 43902 

-- 11. เช็คจำนวน compound_key จาก chembl ที่มีอยู่ใน npass_display_name ด้วย
SELECT COUNT(DISTINCT c.compound_key) 
FROM hdi.hdi_compound_name c 
JOIN hdi.npass_display_name n 
    ON c.compound_key = n.compound_key 
WHERE c.source = 'chembl'; 
-- ผลลัพธ์: 19250 

-- 12. สร้างตารางดึงชื่อ pref_name เฉพาะของ chembl
CREATE TABLE hdi.n_chembl_pref AS 
SELECT DISTINCT ON (compound_key) 
    compound_key, 
    name AS pref_name 
FROM hdi.hdi_compound_name 
WHERE source = 'chembl' AND name_type = 'pref_name' 
ORDER BY compound_key; 

-- 13. สร้างตารางรวมชื่อ display_name (นำข้อมูลจาก np และ chembl มารวมกัน)
CREATE TABLE hdi.n_unified_display_name AS 
SELECT 
    COALESCE(c.compound_key, n.compound_key) AS compound_key, 
    COALESCE(
        c.pref_name, 
        n.pref_name, 
        n.iupac_name, 
        COALESCE(c.compound_key, n.compound_key)
    ) AS display_name, 
    CASE 
        WHEN c.pref_name IS NOT NULL THEN 'chembl_pref' 
        WHEN n.pref_name IS NOT NULL THEN 'np_pref' 
        WHEN n.iupac_name IS NOT NULL THEN 'np_iupac' 
        ELSE 'compound_key_fallback' 
    END AS display_source 
FROM hdi.npass_display_name n 
FULL OUTER JOIN hdi.n_chembl_pref c 
    ON c.compound_key = n.compound_key; 
-- ผลลัพธ์: 225019 

-- 14. เช็คว่ามี compound_key ซ้ำในตารางที่รวมข้อมูลแล้วหรือไม่
SELECT compound_key, COUNT(*) 
FROM hdi.n_unified_display_name 
GROUP BY compound_key 
HAVING COUNT(*) > 1; 
-- ผลลัพธ์: 0 

-- 15. เช็คว่ามีรายการไหนที่ไม่มี display_name หลังจากรวมแล้วหรือไม่
SELECT COUNT(*) 
FROM hdi.n_unified_display_name 
WHERE display_name IS NULL; 
-- ผลลัพธ์: 0 

-- 16. เช็คว่ามีรายการไหนที่ไม่มี compound_key หรือไม่
SELECT COUNT(*) 
FROM hdi.n_unified_display_name 
WHERE compound_key IS NULL; 
-- ผลลัพธ์: 0 

-- 17. เช็คจำนวน compound_key ของ chembl ที่เป็น pref_name
SELECT COUNT(DISTINCT compound_key) 
FROM hdi.hdi_compound_name 
WHERE source = 'chembl' AND name_type = 'pref_name'; 
-- ผลลัพธ์: 42299 

-- 18. เช็คจำนวนรายการที่ใช้ compound_key เป็นชื่อแสดงผลแทนชื่อจริง
SELECT COUNT(*) 
FROM hdi.npass_display_name 
WHERE display_name = compound_key; 
-- ผลลัพธ์: 3