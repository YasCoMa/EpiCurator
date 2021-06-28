# Epitopes Curation Tool

This tool consists of a set of filtering steps to enhance the discovery of new epitopes.

## Requirements:
* External softwares needed:
	- ncbi-blast
	- mafft
	- BepiPred2.0
	
* Python packages needed:
	- os
	- scipy
	- numpy
	- Levenshtein
	- xml
	- statistics

## Usage Instructions
* Preparation:
	1. ````git clone https://github.com/YasCoMa/EpiCurator.git````
	2. ````cd EpiCurator````
	3. ````unzip human_38_proteins````
	4. ````move human_38_proteins.faa into the same directory as epicurator.py````
	5. ````unzip lineages````
	6. ````move lineages folder files into the same directory as epicurator.py````
	7. ````unzip blastdb````
	8. ````Go to EpiMiner folder to see its specific instructions, if you want to use it````

* Pipeline parameters:
	- [function] <br>
		Use to indicate the function you want to execute: <br>
		1 - Analyze mutation details such as positions of AA modification <br>
		2 - Compare overlapping and similarity between epitopes of distinct classes (I, II and B-cell)
		3 - Validate epitopes occurrence in human genome and proteins
		4 - Analyze epitopes provided by pdb structure and Discotope
		5 - Analyze conservation of epitopes in lineages genome sequences and portions of genomes that they were in low, moderated or high coverage
		6 - Check epitopes that are already published in IEDB database
		7 - Calculate correlation between protein and epitope sizes
		8 - Get epitopes coordinates in reference proteins
		9 - Run bepipred
		10 - Analyze output bepipred
		11 - Get epitopes coordinates in genome

	- [file] <br>
		File containing the epitopes separated by comma, with header in the following order (identification, epitope sequence and original protein/epitopes class (for function 2)) <br>
		This parameter is used in all functions except 4, 9 and 10
	
	- [folder] <br>
		Used in functions 4, 9 and 10<br>

* Running modes examples:
	1. Run mode 1: <br>
	````python3 epicurator.py [function] [file] ````

	2. Run mode 2: <br>
	````python3 epicurator.py [function] [folder] ````

## Reference
Article in process.

## Bug Report
Please, use the [Issue](https://github.com/YasCoMa/Epicurator/issues) tab to report any bug.
