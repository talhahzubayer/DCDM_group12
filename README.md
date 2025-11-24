# Data Cleaning and Data Management Coursework - Group 12

## Project Structure

```
data/
├── processed/
│   ├── Disease_information_cleaned.csv         # Disease information metadata file converted to CSV and cleaned
│   ├── GROUP12_clean_table_final.csv           # Final CSV containing all phenotype analysis results data (cleaned)
│   ├── IMPC_parameter_description_cleaned.csv  # IMPC parameter description metadata file converted to CSV and cleaned
│   └── IMPC_procedure_cleaned.csv              # IMPC procedure metadata file converted to CSV and cleaned
│
└── raw/
    ├── phenotype_analysis_results/          # Phenotype analysis result CSV files egressed from TRE
    │   ├── 0000b4e33ek4f65.csv
    │   ├── 000k937al220z48.csv
    │   ├── ...
    │   └── 00b6c53femc2397.csv
    │   (28,590 CSV files in total)
    │
    └── phenotype_analysis_results_transposed/  # Transposed phenotype analysis CSV files for data consolidation
        ├── transposed_0000b4e33ek4f65.csv
        ├── transposed_000k937al220z48.csv
        ├── ...
        └── transposed_00d67j5s7k89yfz.csv
        (28,590 CSV files in total)

database/
├── DCDM_CW1_GROUP12.png                     # Database schema diagram
├── DCDM_CW1_GROUP12.sql                     # Database dump file
├── SQL_create_database.sql                  # Database creation script
└── SQL_query_genes.sql                      # Query gene extraction script

metadata/
├── Disease_information.txt                  # Disease information metadata egressed from TRE
├── IMPC_parameter_description.txt           # IMPC parameter description metadata egressed from TRE
├── IMPC_procedure.txt                       # IMPC procedure metadata egressed from TRE
├── IMPC_SOP.csv                             # IMPC Standard Operating Procedures egressed from TRE
└── query_genes.csv                          # Query genes list egressed from TRE

scripts/
├── 1_transpose_data.R                       # Data transposition script
├── 2_data_combining_and_cleaning.R          # Data combining and cleaning script
├── 3_disease_information_cleaning.R         # Disease information cleaning script
├── 4_parameter_description_cleaning.R       # Parameter description cleaning script
├── 5_procedure_cleaning.R                   # Procedure cleaning script
├── data_loader_module.R                     # Data loader module for impc_dashboard script
└── impc_dashboard.R                         # RShiny dashboard application
```
