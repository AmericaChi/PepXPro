# Predicted and Experimental Peptide Binding Information (PEPBI) Database

---

[Access this dataset on Dryad](https://doi.org/10.5061/dryad.wstqjq2wk)

This dataset is the publicly available data associated with *A Paired Database of Predicted and Experimental Protein Peptide Binding Information*, submitted to *Scientific Data* and hereafter referred to as the PEPBI manuscript. The Predicted and Experimental Peptide Binding Information (PEPBI) dataset provides curated thermodynamic and structural data for protein–peptide complexes. It includes:

* **Experimental Data**: Isothermal titration calorimetry (ITC) data collected through an extensive literature review, compiled in the `PEPBI.xlsx` Excel spreadsheet.
* **Computational Data**: Computationally predicted protein–peptide structures in `.pdb` format, each derived from a shared reference structure (template). These templates are experimentally determined structures from the Protein Data Bank (PDB) and are used as the foundation for generating computational models of mutated variants when experimental structures for those variants are not available. Variants generated from the same template are grouped together and referred to as a binding group. All predicted structures were analyzed using Rosetta’s Interface Analyzer (RIA), and the resulting binding interface metrics are compiled in the PEPBI.xlsx spreadsheet. Structural modifications used to generate the models are documented in accompanying .docx files.
* **Technical Validation Data**: Binding energy distribution data from external protein–peptide and antibody–antigen benchmarks used to evaluate PEPBI predictions. These are compiled in the `Technical Validation.csv` spreadsheet.
* **Utility Script**: A Python script, `PEPBI.py`, that enables users to convert the dataset into a `.csv` file or load it directly as a pandas DataFrame for analysis.

This dataset is intended to support research in protein–peptide interaction prediction, computational modeling, and thermodynamic characterization.

---

## Description of the Data and File Structure

The PEPBI database consists of **A**. an excel spreadsheet containing sequence information, experimental thermodynamic measurements, and computationally predicted binding interaction metrics, **B**. computationally predicted structural models, **C**. an excel spreadsheet containing technical validation data, and **D**. supporting documentation. Below is a description of each component in the repository, along with an overview of how the data is organized.

### 1. PEPBI Spreadsheet (`PEPBI.xlsx`)

This Excel file contains curated data for the protein-peptide complexes included in the PEPBI database. Each row represents a unique complex. To aid interpretation, the contents of the spreadsheet are organized into three general classifications: Sequence Information, Experimental Data, and Computational Data. These classifications serve as descriptive groupings to help users navigate the file structure and are not explicit headers in the spreadsheet itself. This section is organized as follows:

1. First-level bullet points represent the three broad classifications.
2. Second-level bullet points describe logical groupings of related data (e.g., Protein–Peptide Complex Identifiers, Protein Sequence Information, Peptide Sequence Information, etc.), which correspond to the first row of headers within the spreadsheet.
3. Third-level bullet points list the individual column names found under each grouping.
4. Any additional indentation provides supporting information such as formatting details or naming conventions related to column entries.

* **Sequence Information**: Descriptive information for each protein–peptide complex, including identifiers, sequence-level details, and information about its relationship to the reference structure. This metadata helps users interpret the data by specifying what each complex is and how it connects to other entries in the database. The following columns are presented in the order they appear in the spreadsheet:
  * **Protein–Peptide Complex Identifiers**:
    * **Binding Group**: A set of protein–peptide complexes derived from the same template structure, each representing a variant of that reference.
      * **Naming Convention**
        * `Protein - Peptide`, where:
        * `Protein` = protein name
        * `Peptide` = peptide name
        * **Example**: Fyn SH3 - P2L
    * **PEPBI Complex Name**: A unique identifier assigned to each complex in the dataset. These names correspond directly to the filenames of the computationally predicted `.pdb` structures. The naming convention uses the protein and peptide names (often abbreviated) separated by underscores, followed by an index used to distinguish each variant within the same binding group.
      * **Naming Convention**
        * `protein_peptide_#`, where:
        * `protein` = protein name
        * `peptide` = peptide name
        * `#` = index
        * **Example**: fyn_p2l_1
    * **PDB ID**: The accession code for each entry's template structure in the Protein Data Bank (PDB).
    * **Crystallographic Unit**: The number of protein/peptide copies present in the crystal structure. Each crystallographic unit is labeled with `A` for the protein and `B` for the peptide, each followed by the number of copies present in the structure.
      * **Example**: A1B1
    * **Unit Copy**: Indicates which instance of the crystallographic unit the structure represents when multiple versions of the same unit are available. Entries are labeled `C#` for each alternative copy (e.g., `C1`, `C2`), or `X` if no alternative copies exist.
  * **Protein Sequence Information**:
    * **Protein Name**: The name of the protein as reported in the source manuscript.
    * **Protein Sequence**: The protein’s amino acid sequence using one-letter code.
    * **Protein Length**: The total number of amino acids in the protein sequence.
    * **Protein Tag**: Any experimental tags (e.g., `His-tag`, `GST-tag`) that were present on the protein during the ITC experiments.
  * **Peptide Sequence Information**:
    * **Peptide Name**: The name of the peptide as reported in the source manuscript.
    * **Peptide Sequence**: The peptide's amino acid sequence using one-letter code.
    * **Peptide Length**: The total number of amino acids in the peptide sequence.
    * **Peptide Motif**: A recurring peptide sequence important for binding as reported in the source manuscript.
  * **Protein-Peptide Complex Variations**:
    * **Change from Binding Group Reference Complex**: General classifications of how each complex differs from its template structure. These classifications may appear individually or in combination when there are multiple changes needed. If no changes were made, the entry is marked with an `X`. Possible 'change' classifications include:
      * `ITC Temperature (T)`
      * `Peptide Sequence`
      * `Protein Sequence`
      * `Binding Site`
    * **Binding Site**: Indicates which binding site the peptide occupies when the protein has multiple known binding sites. Entries are labeled as `Site 1`, `Site 2`, etc. If only one binding site is known or applicable, the entry is marked as `X`. Documentation explaining the distinction between binding sites is provided in a later section (*see 2. Computational Structures (`Binding Group Structures.zip`)*).
    * **Protein Mutation(s)**: Describes mutations/modifications made to the protein sequence relative to the template structure.
      * **Single-site mutation**: `A5L` (alanine at position 5 mutated to leucine)
      * **Insertion**: `ins940QEPE` (insertion of "QEPE" beginning at position 940)
      * **Single deletion**: `T20del` (threonine at position 20 deleted)
      * **Range deletion**: `F32_A34del` (residues from position 32 to 34 deleted)
      * **Multiple modifications**: Listed and separated by commas (e.g., `A5L, T20del`)
    * **Peptide Mutation(s)**: Describes mutations made to the peptide sequence using the same annotation format as for protein mutations.
* **Experimental Data**: Thermodynamic parameters obtained from isothermal titration calorimetry (ITC) experiments, compiled through an extensive literature review. All values have been converted to consistent units for comparability across entries. The following columns appear in the order presented in the spreadsheet:
  * **Experimental Isothermal Titration Calorimetry (ITC) Data**:
    * **ΔG (kcal/mol)**: Gibbs free energy
    * **ΔG SD (kcal/mol)**: Standard deviation of ΔG
    * **KD (M)**: Dissociation constant
    * **KD SD (M)**: Standard deviation of KD
    * **ΔH (kcal/mol)**: Enthalpy
    * **ΔH SD (kcal/mol)**: Standard deviation of ΔH
    * **TΔS (kcal/mol)**: Temperature-scaled entropy contribution
    * **TΔS SD (kcal/mol)**: Standard deviation of TΔS
    * **ΔS (kcal/mol·K)**: Entropy
    * **ΔS SD (kcal/mol·K)**: Standard deviation of ΔS
    * **N**: Stoichiometry
    * **N SD**: Standard deviation of N
    * **T (K)**: Temperature in Kelvin
  * **Thermodynamic Calculations from ITC Data**: These values were calculated from experimentally reported ITC measurements using standard thermodynamic equations.
    * **Calculated ΔG (kcal/mol)**: Gibbs free energy, calculated as:\
      ΔG = –RTln(1/KD)\
      where R = 0.001987 kcal/mol·K
    * **Calculated TΔS (kcal/mol)**: Temperature-scaled entropy contribution, calculated as:\
      TΔS = ΔH – ΔG
    * **Calculated ΔS (kcal/mol·K)**: Entropy change, calculated as:\
      ΔS = TΔS/T
* **Computational Data**: Includes 40 binding interface metrics generated using Rosetta’s Interface Analyzer (RIA). The following columns are direct outputs from RIA and are presented in the order they appear in the spreadsheet. Descriptions of these metrics can be found in **Table 2: The RIA-calculated Properties** of the corresponding PEPBI manuscript (can also be found in Supplemental Information Zenodo link).
  * **Computational Rosetta Interface Analyzer (RIA) Data**:
    * total_score
    * complex_normalized
    * dG_cross
    * dG_cross/dSAASAx100
    * dG_separated
    * dG_separated/dSASAx100
    * dSASA_hphobic
    * dSASA_int
    * dSASA_polar
    * delta_unsatHbonds
    * dslf_fa13
    * fa_atr
    * fa_dun
    * fa_elec
    * fa_intra_rep
    * fa_intra_sol_xover4
    * fa_rep
    * fa_sol
    * hbond_E_fraction
    * hbond_bb_sc
    * hbond_lr_bb
    * hbond_sc
    * hbond_sr_bb
    * hbonds_int
    * lk_ball_wtd
    * nres_all
    * nres_int
    * omega
    * p_aa_pp
    * packstat
    * per_residue_energy_int
    * pro_close
    * rama_prepro
    * ref
    * sc_value
    * side1_normalized
    * side1_score
    * side2_normalized
    * side2_score
    * yhh_planarity

---

### 2. Computational Structures (`Binding Group Structures.zip`)

This zip folder contains 32 subfolders, each corresponding to a Binding Group in PEPBI as listed in **Table 1: PEPBI Database Summary**. In cases where a binding group contains multiple crystallographic units, additional subfolders are included and labeled according to the corresponding crystallographic unit. Each of these subfolders contains the appropriate .pdb files specific to that unit. Additionally, each binding group subfolder contains a single .docx file that documents the modifications made to generate each computational model from its corresponding template structure. When crystallographic unit subfolders are present, this .docx file remains in the main binding group folder rather than within the individual unit folders.

* **Binding Group Subfolders**: Folder names correspond to the 32 binding groups listed in the ‘Binding Group’ column of the `PEPBI.xlsx` spreadsheet and detailed in *Table 1: PEPBI Database Summary*.
* **Crystallographic Unit Subfolders**: When present, these folders correspond to the ‘Crystallographic Unit’ column of the `PEPBI.xlsx` spreadsheet.
* **Structure Files**: `.pdb` file names match the entries in the ‘PEPBI Complex Name’ column of the `PEPBI.xlsx` spreadsheet.
* **Structure Generation Documents**: `.pdf` file names match those of the 'Binding Group' column of the `PEPBI.xlsx` spreadsheet. Each `.pdf` file includes:
  * A visual representation of the template structure(s)
  * Annotated protein and peptide sequences of the template structure(s), including position numbers
  * A list of all modifications made for each complex (e.g., mutations, insertions, deletions)
  * Highlighted sequences showing the locations of modifications

---

### 3. Technical Validation Data (`Technical Validation.csv`)

The binding energy data used to validate our computational predictions is provided in the `Technical Validation.csv` file. This dataset corresponds to the binding energy distributions shown in **Figure 3**.

This CSV contains a consolidated set of binding energy values from multiple sources:

* **PEPBI Experimental**: Calculated ΔG values (in kcal/mol) from the `Calculated ΔG (kcal/mol)` column in `PEPBI.xlsx`, derived from experimental measurements.
* **PEPBI**: Predicted binding energies (in kcal/mol) using the `dG_separated` metric from Rosetta’s Interface Analyzer (RIA), as provided in the `dG_separated` column of `PEPBI.xlsx`.
* **AB/AG**: Predicted binding energies (in kcal/mol) for mutated antigen–antibody complexes, using the `Mut_IE` metric (equivalent to `dG_separated` from RIA) sourced from [https://doi.org/10.1080/19420862.2024.2440586](https://doi.org/10.1080/19420862.2024.2440586).
* **MM/PBSA**: Binding energies reported as ΔG values (in kcal/mol), converted from the original kJ/mol units, sourced from [https://pubs.acs.org/doi/10.1021/acs.jcim.3c00752](https://pubs.acs.org/doi/10.1021/acs.jcim.3c00752).

> **Note:** This CSV includes only the consolidated binding energy values used for Figure 3. It does not include the accompanying metadata (e.g., PDB identifiers, sequences, mutation details) that were originally reported in the source publications from which the binding energy data were obtained. The data values are listed in the same order as presented in the original sources to preserve consistency.

---

### 4. PEPBI Utility Script (`PEPBI.py`)

`PEPBI.py` is a Python script designed to assist users in accessing and working with the PEPBI dataset by converting it from an Excel `.xlsx` format into a `.csv` file or loading it directly as a pandas DataFrame for further analysis. The script performs the following tasks:

1. **Loading the Spreadsheet into a DataFrame**: The script loads the `PEPBI.xlsx` file into a pandas DataFrame, skipping the first row (which contains merged group labels) and using the second row for column headers.
2. **Data Preview and Information**: After loading the data into the DataFrame, the script previews the first 10 rows and displays the number of rows and columns in the dataset.
3. **Symbol Replacement**: It replaces certain symbols (Δ, α, and β) in the column names and data with their corresponding letters (d, a, and b) for consistency and compatibility with CSV formatting.
4. **Saving to CSV**: Finally, the modified DataFrame is saved as a `.csv` file.

The script only requires the `pandas` library and provides a simple way to prepare and save the dataset for analysis. It helps users easily access and manipulate the PEPBI dataset.
