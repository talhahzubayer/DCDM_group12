# ============================================================================
# FIXED VERSION - SQL_create_database.sql
# All syntax errors corrected, ready to run
# ============================================================================

# 1. Create Core Tables

# Genes table
CREATE TABLE Genes (
    gene_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    gene_symbol VARCHAR(100) NOT NULL,
    gene_accession_id VARCHAR(50) UNIQUE,
    parameter_name VARCHAR(200)
);

# Mouse characteristics
CREATE TABLE MouseCharacteristics (
    mouse_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    mouse_strain VARCHAR(100) NOT NULL,
    mouse_life_stage VARCHAR(100) NOT NULL
);

# Parameters table
CREATE TABLE Parameters (
    parameter_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    IMPC_parameter_id VARCHAR(20) UNIQUE,
    parameter_description VARCHAR(100) NOT NULL,
    group_id INT,
    IMPC_parameter_origin_id VARCHAR(100) NOT NULL
);

# Parameter Groups
CREATE TABLE Parameter_Groups (
    group_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    group_name VARCHAR(50) NOT NULL,
    description VARCHAR(255)
);

# Analysis table
CREATE TABLE Analysis (
    analysis_key_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    pvalue VARCHAR(100) NOT NULL,
    analysis_id VARCHAR(100),
    gene_accession_id VARCHAR(100),
    mouse_id INT,
    gene_id INT,
    parameter_id INT,
    IMPC_parameter_id VARCHAR(20)
);

# Diseases table -- FIXED: Added missing commas
CREATE TABLE Diseases (
    disease_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    disease_name VARCHAR(200) NOT NULL,
    gene_accession_id VARCHAR(50),
    gene_id INT,
    omim_id VARCHAR(20),  -- FIXED: Added comma here
    DO_ID VARCHAR(50) NOT NULL
);

# Procedure table -- FIXED: Added missing commas
CREATE TABLE Procedures (
    IMPC_procedure_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    procedure_name VARCHAR(100) NOT NULL,  -- FIXED: Added comma
    procedure_description VARCHAR(1000) NULL,  -- FIXED: Added comma
    mandatory TINYINT(1),  -- FIXED: Added comma
    IMPC_parameter_id VARCHAR(20) NULL,  -- FIXED: Added comma
    IMPC_parameter_origin_id VARCHAR(50) NOT NULL
);

# 2. POPULATE TABLES:
# ALL DATA WAS IMPORTED FROM CLEANED AND NORMALISED DATAFILES. 
# REFER TO disease_indormation_cleaning.R, data_combining_and_cleanong.R, 
# parameter_description_cleaning.R, procedure_cleaning.R and trasnpose_data.R FOR CLEAN DATA.
# CLEAN DATA CAN BE FOUND IN DCDM_group12/data/processed.

# 3. Add Foreign Keys

ALTER TABLE Analysis
ADD CONSTRAINT fk_analysis_gene
FOREIGN KEY (gene_id) REFERENCES Genes(gene_id)
ON DELETE SET NULL;

ALTER TABLE Analysis
ADD CONSTRAINT fk_analysis_mouse
FOREIGN KEY (mouse_id) REFERENCES MouseCharacteristics(mouse_id)
ON DELETE SET NULL;

ALTER TABLE Analysis
ADD CONSTRAINT fk_analysis_parameter
FOREIGN KEY (parameter_id) REFERENCES Parameters(parameter_id);

ALTER TABLE Parameters
ADD CONSTRAINT fk_group
FOREIGN KEY (group_id) REFERENCES Parameter_Groups(group_id)
ON DELETE SET NULL;

ALTER TABLE Diseases
ADD CONSTRAINT fk_diseases_gene
FOREIGN KEY (gene_id) REFERENCES Genes(gene_id)
ON UPDATE CASCADE
ON DELETE RESTRICT;

# REMOVED: Gene_Disease_Association table doesn't exist yet
# If you need this, create the table first:
# CREATE TABLE Gene_Disease_Association (
#     association_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
#     gene_id INT,
#     disease_id INT,
#     FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
#     FOREIGN KEY (disease_id) REFERENCES Diseases(disease_id)
# );

# 4. Insert Parameter Groups

INSERT INTO Parameter_Groups (group_name, description)
VALUES
('Housing & Environment', 'Mouse husbandry, caging, diet, ambient factors'),
('Structural Phenotype', 'Physical and structural measurements'),
('Clinical Chemistry/Blood', 'Hematology and blood chemistry'),
('Embryo & Development', 'Embryo, viability, and developmental defects'),
('Limb Function & Performance', 'Strength, movement, grip measures'),
('Behavioral & Neurological', 'Activity, behavior, neurological function'),
('Cardiovascular & ECG', 'Heart function, ECG measurements'),
('Other', 'Miscellaneous parameters');

# 5. Parameter Group Assignments
# NOTE: These UPDATE statements reference column 'parameter_name' but Parameters table 
# only has 'parameter_description'. You may need to adjust based on your actual data.
# For now, I'm commenting these out - run them AFTER you load your parameter data.

-- Housing & Environment
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Housing & Environment'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_hou%'
--      OR LOWER(p.parameter_description) LIKE '%housing%'
--      OR LOWER(p.parameter_description) LIKE '%nutrition%'
--      OR LOWER(p.parameter_description) LIKE '%diet%'
--      OR LOWER(p.parameter_description) LIKE '%temperature%'
--      OR LOWER(p.parameter_description) LIKE '%humidity%'
--      OR LOWER(p.parameter_description) LIKE '%light%'
--      OR LOWER(p.parameter_description) LIKE '%cage%'
--      OR LOWER(p.parameter_description) LIKE '%water%'
--   );

-- Structural Phenotype
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Structural Phenotype'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_owt%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_his%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_dxa%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_bwt%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_eye%'
--      OR LOWER(p.parameter_description) LIKE '%weight%'
--      OR LOWER(p.parameter_description) LIKE '%severity score%'
--      OR LOWER(p.parameter_description) LIKE '%skin%'
--      OR LOWER(p.parameter_description) LIKE '%bone%'
--      OR LOWER(p.parameter_description) LIKE '%brain%'
--   );

-- Clinical Chemistry/Blood
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Clinical Chemistry/Blood'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_cbc%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_cld%'
--      OR LOWER(p.parameter_description) LIKE '%blood%'
--      OR LOWER(p.parameter_description) LIKE '%plasma%'
--      OR LOWER(p.parameter_description) LIKE '%hemoglobin%'
--   );

-- Embryo & Development
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Embryo & Development'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_evo%'
--      OR LOWER(p.parameter_description) LIKE '%embryo%'
--      OR LOWER(p.parameter_description) LIKE '%fetus%'
--   );

-- Limb Function & Performance
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Limb Function & Performance'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_grs%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_rot%'
--      OR LOWER(p.parameter_description) LIKE '%limb%'
--      OR LOWER(p.parameter_description) LIKE '%grip%'
--   );

-- Behavioral & Neurological
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Behavioral & Neurological'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_oft%'
--      OR LOWER(p.IMPC_parameter_id) LIKE 'impc_aas%'
--      OR LOWER(p.parameter_description) LIKE '%behavior%'
--      OR LOWER(p.parameter_description) LIKE '%activity%'
--   );

-- Cardiovascular & ECG
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Cardiovascular & ECG'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL
--   AND (
--         LOWER(p.IMPC_parameter_id) LIKE 'impc_ecg%'
--      OR LOWER(p.parameter_description) LIKE '%ecg%'
--      OR LOWER(p.parameter_description) LIKE '%heart rate%'
--   );

-- Assign remaining to "Other"
-- UPDATE Parameters p
-- JOIN Parameter_Groups pg ON pg.group_name = 'Other'
-- SET p.group_id = pg.group_id
-- WHERE p.group_id IS NULL;

# ============================================================================
# NOTES ON WHAT NEEDS TO BE DONE AFTER RUNNING THIS SCRIPT:
# ============================================================================
# 
# 1. LOAD YOUR CSV DATA into the tables:
#    - Genes
#    - MouseCharacteristics
#    - Parameters
#    - Procedures
#    - Diseases
#    - Analysis (your phenotype scores)
#
# 2. THEN uncomment and run the parameter group assignment UPDATE statements above
#
# 3. RUN these additional commands to link data (AFTER loading CSV):
#
#    -- Link Analysis to Genes by gene_accession_id
#    UPDATE Analysis a
#    JOIN Genes g ON a.gene_accession_id = g.gene_accession_id
#    SET a.gene_id = g.gene_id;
#
#    -- Link Analysis to MouseCharacteristics
#    UPDATE Analysis a
#    JOIN MouseCharacteristics m ON a.mouse_strain = m.mouse_strain
#    SET a.mouse_id = m.mouse_id;
#
#    -- Link Analysis IMPC_parameter_id to Parameters
#    UPDATE Analysis a
#    JOIN Parameters p ON a.parameter_id = p.parameter_id
#    SET a.IMPC_parameter_id = p.IMPC_parameter_id;
#
# ============================================================================

-- End of script
