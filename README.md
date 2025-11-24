# Data Cleaning and Data Management Coursework - Group 12

## Project Structure

```
data/
├── processed/                                  # CSV files that have been generated after running our R scripts
│   ├── Disease_information_cleaned.csv         # Disease information metadata file converted to CSV and cleaned
│   ├── GROUP12_clean_table_final.csv           # Final CSV containing all phenotype analysis results data (cleaned)
│   ├── IMPC_parameter_description_cleaned.csv  # IMPC parameter description metadata file converted to CSV and cleaned
│   └── IMPC_procedure_cleaned.csv              # IMPC procedure metadata file converted to CSV and cleaned
│
└── raw/
    ├── phenotype_analysis_results/             # Phenotype analysis result CSV files egressed from TRE
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

## Setting Up Local Database with DBeaver

### Prerequisites
- **DBeaver Community Edition**
- **MySQL Server** installed and running on your local machine
- The `DCDM_CW1_GROUP12.sql` database dump file from the `database/` folder

### Step-by-Step Instructions

#### Step 1: Install and Launch DBeaver
1. Download and install DBeaver from https://dbeaver.io/download/](https://dbeaver.io/download/
2. Launch DBeaver

#### Step 2: Create a MySQL Connection
1. Click on **Database** → **New Database Connection** (or click the plug icon)
2. Select **MySQL** from the list of databases
3. Click **Next**
4. Enter your MySQL connection details:
   - **Host**: `localhost`
   - **Port**: `3306` (default MySQL port)
   - **Database**: Leave empty for now
   - **Username**: Your MySQL username (e.g., `root`)
   - **Password**: Your MySQL password
5. Click **Test Connection** to verify the connection works
   - If prompted, download the MySQL driver
6. Click **Finish**

#### Step 3: Create an Empty Database
**IMPORTANT**: You must create an empty database before importing the dump file.

1. In the **Database Navigator** panel (left side), expand your MySQL connection
2. Right-click on **Databases** → **Create New Database**
3. Enter the database name: `DCDM_CW1_GROUP12`
4. Set character set to: `utf8mb4`
5. Set collation to: `utf8mb4_general_ci`
6. Click **OK**

#### Step 4: Import the Database Dump File
1. Right-click on the newly created `DCDM_CW1_GROUP12` database
2. Select **Tools** → **Restore Database**
3. In the restore dialog, click **Browse** to locate the dump file
4. Navigate to your project's `database/` folder and select `DCDM_CW1_GROUP12.sql`
5. Click **Start** to begin the restoration process
6. Wait for the import to complete (this may take a few minutes)
7. You should see a success message indicating the database was restored successfully

#### Step 5: Verify the Import
1. In the **Database Navigator**, expand the `DCDM_CW1_GROUP12` database
2. Expand the **Tables** folder
3. You should see all 6 tables:
   - `Analysis`
   - `Genes`
   - `Parameters`
   - `Parameter_Groups`
   - `Procedures`
   - `Diseases`
4. Right-click on any table and select **View Data** to verify records were imported

---

## Connecting Local Database to RShiny Dashboard

### Prerequisites
- Local MySQL database set up with DBeaver (see previous section)
- R and RStudio installed
- Required R packages:
  - `tidyverse`       - includes dplyr, tidyr, readr, etc.
  - `RMySQL`          - for database connectivity
  - `shiny`           - for dashboard interface
  - `shinydashboard`  - for dashboard layout
  - `plotly`          - for interactive plots
  - `DT`              - for interactive tables
  - `umap`            - for gene clustering
  - `metap`           - for Fisher's method p-value combination

You can install all required packages at once with:
```r
install.packages(c("tidyverse", "RMySQL", "shiny", "shinydashboard", 
                   "plotly", "DT", "umap", "metap"))
```

### Step-by-Step Instructions

#### Step 1: Get Your Database Credentials from DBeaver
1. In DBeaver, right-click on your MySQL connection
2. Select **Edit Connection**
3. Note down the following information:
   - **Host**: (typically `localhost` or `127.0.0.1`)
   - **Port**: (typically `3306`)
   - **Database**: `DCDM_CW1_GROUP12`
   - **Username**: Your MySQL username
   - **Password**: Your MySQL password

#### Step 2: Update the Data Loader Module
1. Open the `scripts/data_loader_module.R` file in RStudio or any text editor
2. Locate the configuration section at the top of the file
3. Find the following variables and update them with your credentials:

```r
# CONFIGURATION 
# Data source: "csv" or "database"
DATA_SOURCE <- "database" # Make sure it is set to database

# Database configuration
DB_CONFIG <- list(
  host = "localhost",              # Change if your MySQL server is on a different host
  port = 3306,					   # Default MySQL port
  dbname = "DCDM_CW1_GROUP12",     # Database name should match what you created in DBeaver
  user = "your_username",          # Replace with your MySQL username from DBeaver
  password = "your_password"       # Replace with your MySQL password from DBeaver
)
```

#### Step 3: Save and Test the Connection
1. Save the `data_loader_module.R` file
2. Open `scripts/impc_dashboard.R` in RStudio
3. Run the dashboard:
   ```r
   shiny::runApp("scripts/impc_dashboard.R")
   ```
4. If the connection is successful, the dashboard should load data from your local MySQL database
5. If you encounter connection errors, verify:
   - MySQL server is running
   - Database credentials are correct
   - The `DCDM_CW1_GROUP12` database exists and contains data
   - Required R packages (`RMySQL`, `DBI`) are installed

#### Step 5: Switching Between CSV and Database Mode
To switch back to using CSV files instead of the database:
1. Open `scripts/data_loader_module.R`
2. Change `DATA_SOURCE <- "database"` to `DATA_SOURCE <- "csv"`
3. Save the file and restart the dashboard


## RShiny Dashboard

### Gene Analysis Tab
<img width="1919" height="804" alt="image" src="https://github.com/user-attachments/assets/4319c021-293a-46e1-8d2f-c72726240a49" />

### Phenotype Analysis Tab
<img width="1919" height="833" alt="image" src="https://github.com/user-attachments/assets/2e289b44-ca3d-4c07-ab47-bcf3110bd9d8" />

### Gene Clustering Tab
<img width="1919" height="928" alt="image" src="https://github.com/user-attachments/assets/09e1a960-90fa-4687-b6e2-2fe9bca68745" />

### Four Query Genes Tab
<img width="1919" height="973" alt="image" src="https://github.com/user-attachments/assets/800a7d86-fbae-40fb-9e1b-619979e85a64" />







