
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
# Diseases table
CREATE TABLE Diseases (
    disease_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    disease_name VARCHAR(200) NOT NULL,
    gene_accession_id VARCHAR(50),
    gene_id INT,
    omim_id VARCHAR(20)
    DO_ID VARCHAR(50) NOT NULL
);

# Procedure table
CREATE TABLE Procedures (
    IMPC_procedure_id INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
    procedure_name VARCHAR(100) NOT NULL
    procedure_description VARCHAR(1000) NULL
    mandatory TINYINT(1)
    IMPC_parameter_id VARCHAR(20) NULL
    IMPC_parameter_origin_id VARCHAR(50) NOT NULL
);
# 2. POPULATE TABLES:
# ALL DATA WAS IMPORTED FROM CLEANED AND NORMALISED DATAFILES. 
# REFER TO disease_indormation_cleaning.R, data_combining_and_cleanong.R, parameter_description_cleaning.R, procedure_cleaning.R and trasnpose_data.R FOR CLEAN DATA.
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

INSERT INTO Analysis (pvalue, analysis_id, gene_accession_id)
SELECT pvalue, analysis_id, gene_accession_id
FROM Genes;


-- Link Analysis to Genes by gene_accession_id
UPDATE Analysis a
JOIN Genes g ON a.gene_accession_id = g.gene_accession_id
SET a.gene_id = g.gene_id;

-- Link Analysis to MouseCharacteristics
UPDATE Analysis a
JOIN MouseCharacteristics m ON a.mouse_strain = m.mouse_strain
SET a.mouse_id = m.mouse_id;


ALTER TABLE Parameters
ADD CONSTRAINT fk_group
FOREIGN KEY (group_id) REFERENCES Parameter_Groups(group_id)
ON DELETE SET NULL;

ALTER TABLE Diseases
ADD CONSTRAINT fk_diseases_gene
FOREIGN KEY (gene_id) REFERENCES Genes(gene_id)
ON UPDATE CASCADE
ON DELETE RESTRICT;

ALTER TABLE Gene_Disease_Association
ADD CONSTRAINT fk_gene
FOREIGN KEY (gene_id) REFERENCES Genes(gene_id),
ADD CONSTRAINT fk_disease
FOREIGN KEY (disease_id) REFERENCES Diseases(disease_id);


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

-- Housing & Environment
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Housing & Environment'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_hou%'
     OR LOWER(p.parameter_name) LIKE '%housing%'
     OR LOWER(p.parameter_name) LIKE '%nutrition%'
     OR LOWER(p.parameter_name) LIKE '%diet%'
     OR LOWER(p.parameter_name) LIKE '%temperature%'
     OR LOWER(p.parameter_name) LIKE '%humidity%'
     OR LOWER(p.parameter_name) LIKE '%light%'
     OR LOWER(p.parameter_name) LIKE '%cage%'
     OR LOWER(p.parameter_name) LIKE '%water%'
  );




-- Structural Phenotype
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Structural Phenotype'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_owt%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_his%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_dxa%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_bwt%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_eye%'
     OR LOWER(p.parameter_name) LIKE '%weight%'
     OR LOWER(p.parameter_name) LIKE '%severity score%'
     OR LOWER(p.parameter_name) LIKE '%skin%'
     OR LOWER(p.parameter_name) LIKE '%coat%'
     OR LOWER(p.parameter_name) LIKE '%bone%'
     OR LOWER(p.parameter_name) LIKE '%skeleton%'
     OR LOWER(p.parameter_name) LIKE '%muscle%'
     OR LOWER(p.parameter_name) LIKE '%brain%'
     OR LOWER(p.parameter_name) LIKE '%eye%'
     OR LOWER(p.parameter_name) LIKE '%heart%'
     OR LOWER(p.parameter_name) LIKE '%lung%'
     OR LOWER(p.parameter_name) LIKE '%liver%'
     OR LOWER(p.parameter_name) LIKE '%spleen%'
     OR LOWER(p.parameter_name) LIKE '%kidney%'
     OR LOWER(p.parameter_name) LIKE '%stomach%'
     OR LOWER(p.parameter_name) LIKE '%intestine%'
     OR LOWER(p.parameter_name) LIKE '%pancreas%'
     OR LOWER(p.parameter_name) LIKE '%thyroid%'
     OR LOWER(p.parameter_name) LIKE '%adrenal%'
     OR LOWER(p.parameter_name) LIKE '%pituitary%'
     OR LOWER(p.parameter_name) LIKE '%thymus%'
     OR LOWER(p.parameter_name) LIKE '%gland%'
     OR LOWER(p.parameter_name) LIKE '%testis%'
     OR LOWER(p.parameter_name) LIKE '%ovary%'
     OR LOWER(p.parameter_name) LIKE '%uterus%'
     OR LOWER(p.parameter_name) LIKE '%bladder%'
     OR LOWER(p.parameter_name) LIKE '%epididymis%'
     OR LOWER(p.parameter_name) LIKE '%prostate%'
     OR LOWER(p.parameter_name) LIKE '%seminal vesicle%'
     OR LOWER(p.parameter_name) LIKE '%penis%'
     OR LOWER(p.parameter_name) LIKE '%dental%'
     OR LOWER(p.parameter_name) LIKE '%tooth%'
     OR LOWER(p.parameter_name) LIKE '%anatomical%'
     OR LOWER(p.parameter_name) LIKE '%structure%'
     OR LOWER(p.parameter_name) LIKE '%length%'
     OR LOWER(p.parameter_name) LIKE '%thickness%'
     OR LOWER(p.parameter_name) LIKE '%size%'
     OR LOWER(p.parameter_name) LIKE '%dimension%'
     OR LOWER(p.parameter_name) LIKE '%morphology%'
     OR LOWER(p.parameter_name) LIKE '%composition%'
  );


-- Clinical Chemistry/Blood
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Clinical Chemistry/Blood'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_cbc%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_cld%'
     OR LOWER(p.parameter_name) LIKE '%blood%'
     OR LOWER(p.parameter_name) LIKE '%plasma%'
     OR LOWER(p.parameter_name) LIKE '%serum%'
     OR LOWER(p.parameter_name) LIKE '%whole blood%'
     OR LOWER(p.parameter_name) LIKE '%hemoglobin%'
     OR LOWER(p.parameter_name) LIKE '%lymphocyte%'
     OR LOWER(p.parameter_name) LIKE '%neutrophil%'
     OR LOWER(p.parameter_name) LIKE '%erythrocyte%'
     OR LOWER(p.parameter_name) LIKE '%platelet%'
     OR LOWER(p.parameter_name) LIKE '%leukocyte%'
     OR LOWER(p.parameter_name) LIKE '%glucose%'
     OR LOWER(p.parameter_name) LIKE '%cholesterol%'
     OR LOWER(p.parameter_name) LIKE '%protein%'
     OR LOWER(p.parameter_name) LIKE '%acid%'
     OR LOWER(p.parameter_name) LIKE '%hormone%'
     OR LOWER(p.parameter_name) LIKE '%immunoglobulin%'
     OR LOWER(p.parameter_name) LIKE '%assay%'
     OR LOWER(p.parameter_name) LIKE '%bilirubin%'
     OR LOWER(p.parameter_name) LIKE '%alt%'
     OR LOWER(p.parameter_name) LIKE '%enzyme%'
     OR LOWER(p.parameter_name) LIKE '%chemistry%'
     OR LOWER(p.parameter_name) LIKE '%metabolite%'
  );


-- Embryo & Development
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Embryo & Development'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_evo%'
     OR LOWER(p.parameter_name) LIKE '%embryo%'
     OR LOWER(p.parameter_name) LIKE '%placenta%'
     OR LOWER(p.parameter_name) LIKE '%fetus%'
     OR LOWER(p.parameter_name) LIKE '%development%'
     OR LOWER(p.parameter_name) LIKE '%viability%'
     OR LOWER(p.parameter_name) LIKE '%dead%'
     OR LOWER(p.parameter_name) LIKE '%stillborn%'
     OR LOWER(p.parameter_name) LIKE '%gestation%'
     OR LOWER(p.parameter_name) LIKE '%umbilic%'
     OR LOWER(p.parameter_name) LIKE '%yolk sac%'
     OR LOWER(p.parameter_name) LIKE '%defect%'
     OR LOWER(p.parameter_name) LIKE '%anomaly%'
  );

-- Limb Function & Performance
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Limb Function & Performance'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_grs%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_rot%'
     OR LOWER(p.parameter_name) LIKE '%limb%'
     OR LOWER(p.parameter_name) LIKE '%forelimb%'
     OR LOWER(p.parameter_name) LIKE '%hindlimb%'
     OR LOWER(p.parameter_name) LIKE '%grip%'
     OR LOWER(p.parameter_name) LIKE '%paw%'
     OR LOWER(p.parameter_name) LIKE '%digit%'
  );

-- Behavioral & Neurological
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Behavioral & Neurological'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_oft%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_aas%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_ecs%'
     OR LOWER(p.IMPC_parameter_id) LIKE 'impc_het%'
     OR LOWER(p.parameter_name) LIKE '%behavior%'
     OR LOWER(p.parameter_name) LIKE '%activity%'
     OR LOWER(p.parameter_name) LIKE '%locomotor%'
     OR LOWER(p.parameter_name) LIKE '%anxiety%'
     OR LOWER(p.parameter_name) LIKE '%sleep%'
     OR LOWER(p.parameter_name) LIKE '%startle%'
     OR LOWER(p.parameter_name) LIKE '%memory%'
     OR LOWER(p.parameter_name) LIKE '%neurological%'
     OR LOWER(p.parameter_name) LIKE '%brainwave%'
     OR LOWER(p.parameter_name) LIKE '%eeg%'
  );

-- Cardiovascular & ECG
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Cardiovascular & ECG'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL
  AND (
        LOWER(p.IMPC_parameter_id) LIKE 'impc_ecg%'
     OR LOWER(p.parameter_name) LIKE '%ecg%'
     OR LOWER(p.parameter_name) LIKE '%heart rate%'
     OR LOWER(p.parameter_name) LIKE '%cardiac%'
     OR LOWER(p.parameter_name) LIKE '%electrocardiogram%'
     OR LOWER(p.parameter_name) LIKE '%blood pressure%'
  );

-- Step 2: Assign remaining unmatched rows to "Other"
UPDATE DCDM_CW1_GROUP12.Parameters p
JOIN DCDM_CW1_GROUP12.Parameter_Groups pg
  ON pg.group_name = 'Other'
SET p.group_id = pg.group_id
WHERE p.group_id IS NULL;



-- Assuming parameter_id is already populated in Analysis:
UPDATE Analysis a
JOIN Parameters p ON a.parameter_id = p.parameter_id
SET a.IMPC_parameter_id = p.IMPC_parameter_id;

-- Populate MouseCharacteristics from Genes

INSERT INTO MouseCharacteristics (mouse_strain, mouse_life_stage)
SELECT DISTINCT mouse_strain, mouse_life_stage
FROM Genes;



I
