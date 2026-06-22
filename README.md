# miRNA node (Python)

This script relates the differentially expressed miRNAs in the studied datasets to their biological nodes. To do so, it generates a list of all differentially expressed miRNAs identified in each study, along with other relevant information associated with each of them. It then incorporates, when available, the nomenclature equivalences of each mature miRNA across the different databases and it also adds information about the evolutionary node to which each miRNA belongs. Additionally, it generates a table with the differentially expressed miRNAs that have been detected with MirGeneDB, but not with the other databases.

---

## Input files

- **Directory containing sRNAde results (differential expression matrices)**
  - Selected by the user 
  - The folder for the assignment method must be named: 
	-	de_rcsa 
	-	de_rcadj
  - The matrices must be named: 
	-	DESeq2_[Condition 1]_vs_[Condition 2].tsv
	-	edgeR_[Condition 1]_vs_[Condition 2].tsv
	-	limma_[Condition 1]_vs_[Condition 2].tsv

- **Directory containing db_relations file**
  - Selected by the user 

---

## Requirements

- **Python 3**
- Required libraries:
  - `pandas`
  - `matplotlib`
  - `seaborn`

No additional external dependencies are required.

---

## Script usage

Run the script from the terminal or a Python environment:

```bash
python mirna_node.py -i <input_file> -o <output_file> [-s <study>] [-d <db_relations>] [-sp <specie>]
````

Arguments
| Argument | Alternative_Argument | Required | Description | Default |
| --- | --- | --- | --- | --- |
| -i | --input | Yes | Path to the input file | - |
| -o | --output | Yes | Path to the output file | - |
| -s | --study | No | Select the name of the study | SRP |
| -d | --db_relations | No | Path to db_relations | /shared/bak/TFG/serrano/scripts_outputs/db_relations/db_relations.tsv |
| -sp | --specie | No | Specie | hsa |

---

## Output

The script generates five output files:

1. **Tab-separated values file (tsv) containing all the miRNA differentially expressed and their node**

2. **TSV summarizing the number of miRNAs with assigned nodes.**

3. **TSV containing the miRNA differentially expressed only detected with MirGeneDB**

4. **TSV summarizing the number each miRNA is only detected with MirGeneDB**

5. **Hybrid plot combining a box plot and a strip plot that summarizes the percentage of miRNAs with biological nodes (png)**

---

## How the script works

The script is divided into three main stages:

### 1️. File selection and loading into a `DataFrame`

- Reading the selected files into a `DataFrame`.

### 2. Generation of the `DataFrame` containing the relationships between miRNAs and their biological nodes, and generation of the plot.

- Using the differentially detected miRNAs in the studies and db_relations, each miRNA is annotated with its biological origin node, if available, as well as its alternative names across the different databases.
- Store these relationships in a new `DataFrame` and save it as a tsv file.
- Summarize this information for each study and methodological combination by counting the number of differentially expressed miRNAs detected and the number of those miRNAs with an assigned biological node.
- Store the information in a new `DataFrame` and save it as a tsv file.
- Using the last DataFrame, generate a graphical representation that summarizes all the information.

### 3. Generation of the `DataFrame` containing only the miRNAs detected using MirGeneDB.

- Check in the previously generated `DataFrame` which differentially expressed miRNAs are detected using MirGeneDB and not detected using the other databases.
- Store the information in a new `DataFrame` and save it as a tsv file.
- Summarize the information for each database in a new `DataFrame` and save it as a tsv file.
