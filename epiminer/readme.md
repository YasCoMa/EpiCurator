# EpiMiner

The EpiMiner starts getting all papers that mention the epitope and extrating their content using the NCBI [API](https://www.ncbi.nlm.nih.gov/home/develop/api/). These XML files are cleaned removing hypertext markup and references to figures, tables and supplementary materials. The paragraphs of the remaining articles content are processed by Natural language processing steps to extract sentences, tokens, stopwords removal to remove words extremely common in english language and do not help to identify the context of interest, prioritizing tokens using part-of-speech tagging to keep just nouns and verbs. Then, the sentences filtered goes to the task that identifies the epitope sequences in evaluation among the tokens. Finally, a report is made by epitope with the article identifiers and the sentences.

## Requirements:
* Python packages needed:
	- pandas
	- urllib
	- os
	- rdflib
	- mlxtend
	- re
	- string
	- unicodedata
	- inflect
	- nltk
	- Bio python
	- xml
	- itertools
	- lxml
	- bs4 (beautiful soup)

## Usage Instructions

* Pipeline parameters:
	- __-fo__ or __--folder__ <br>
		Folder to store the files (use the folder where the other required file can be found)

	- __-rt__ or __--running_type__ <br>
		Use to indicate which execution step you want to run for mode 1 (it is desirable following the order showed): <br>
		0 (default) - Run all steps <br>
		1 - Run step 1 (Get mentions of epitopes in PMC articles) <br>
		2 - Run step 2 (Get the PMC or Pubmed files, clean and store them) <br>
		3 - Run step 3 (Get the exact sentences where the epitopes were found)

	- __-fp__ or __--file_epitopes__ <br>
		File with the epitopes (one column with the epitope sequence in tsv format, without header, located into the folder assigned in the above parameter)<br>

	- __-fe__ or __--file_evaluation__ <br>
		File exported after step 1 execution in tsv format<br>
		
* Running modes examples:
	1. Running all three steps: <br>
		````python3 epiminer.py -rt 0 -fo data/ -fp epitopes.tsv````
	1. Running from step 2: <br>
		````python3 epiminer.py -rt 2 -fo data/ -fp literature_evaluation_pairs.tsv````

## Reference

## Bug Report
Please, use the [Issue](https://github.com/YasCoMa/Epicurator/issues) tab to report any bug.
