-- Query 1: Basic Retrieval of Gene Details
-- Purpose: Get all columns (gene_id, gene_symbol, accession_id, etc.) 
--          for the gene 'Smarcd3' from the Genes table.
SELECT 
    *
FROM 
    Genes
WHERE 
    gene_symbol = 'Smarcd3';

-- Query 2: Retrieve Gene Details and Associated Analysis Results
-- Purpose: Join the Genes table with the Analysis table to see all 
--          experimental/statistical data (like pvalue and analysis_key_id)
--          associated with 'Smarcd3'.
SELECT
    G.gene_symbol,
    G.gene_accession_id,
    A.analysis_key_id,
    A.pvalue AS analysis_pvalue,
    A.parameter_id
FROM
    Genes G
JOIN
    Analysis A ON G.gene_id = A.gene_id
WHERE
    G.gene_symbol = 'Smarcd3';
# Example Queries

-- Query gene and its analysis results
SELECT G.gene_symbol,
       G.gene_accession_id,
       A.analysis_key_id,
       A.pvalue AS analysis_pvalue,
       P.parameter_description
FROM Genes G
JOIN Analysis A ON G.gene_id = A.gene_id
JOIN Parameters P ON A.parameter_id = P.parameter_id
WHERE G.gene_symbol = 'Smarcd3';

-- Left join to include all genes even if no analysis exists
SELECT G.gene_symbol,
       G.gene_accession_id,
       A.analysis_key_id,
       A.pvalue AS analysis_pvalue,
       P.parameter_description
FROM Genes G
LEFT JOIN Analysis A ON G.gene_id = A.gene_id
LEFT JOIN Parameters P ON A.parameter_id = P.parameter_id
WHERE G.gene_symbol = 'Smarcd3';


UPDATED:
# This query looks up the gene Smarcd3 in the Genes table and then pulls in any related analysis results and parameter descriptions by 
matching gene_id and parameter_id. It uses LEFT JOINs, so the Smarcd3 gene will still appear even if no matching 
analysis or parameter data exists—those fields will simply show up as NULL.

SELECT G.gene_symbol,
       G.gene_accession_id,
       A.analysis_key_id,
       A.pvalue AS analysis_pvalue,
       P.parameter_description
FROM Genes G
LEFT JOIN Analysis A ON G.gene_id = A.gene_id
LEFT JOIN Parameters P ON A.parameter_id = P.parameter_id
WHERE G.gene_symbol = 'Smarcd3';

SELECT G.gene_symbol,
       G.gene_accession_id,
       A.analysis_key_id,
       A.pvalue AS analysis_pvalue,
       P.parameter_description
FROM Genes G
LEFT JOIN Analysis A ON G.gene_id = A.gene_id
LEFT JOIN Parameters P ON A.parameter_id = P.parameter_id
WHERE G.gene_symbol = 'Ppp3cc';

SELECT G.gene_symbol,
       G.gene_accession_id,
       A.analysis_key_id,
       A.pvalue AS analysis_pvalue,
       P.parameter_description
FROM Genes G
LEFT JOIN Analysis A ON G.gene_id = A.gene_id
LEFT JOIN Parameters P ON A.parameter_id = P.parameter_id
WHERE G.gene_symbol = 'Rab12';


SELECT G.gene_symbol,
       G.gene_accession_id,
       A.analysis_key_id,
       A.pvalue AS analysis_pvalue,
       P.parameter_description
FROM Genes G
LEFT JOIN Analysis A ON G.gene_id = A.gene_id
LEFT JOIN Parameters P ON A.parameter_id = P.parameter_id
WHERE G.gene_symbol = 'Klhl33';


SELECT *
FROM Analysis
WHERE gene_id IN (SELECT gene_id FROM Genes WHERE gene_symbol = 'Klhl33');


