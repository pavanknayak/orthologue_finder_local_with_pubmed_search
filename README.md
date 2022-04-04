# Orthologue finder and NCBI searcher for ENSEMBL ID Gene Lists

This python program is an improvement off of the initial Gene list NCBI searcher I created, shown here: [EE282 Project](https://github.com/pnayak93/ee282/blob/main/Final_Project_Writeup.md)

The original project consisted of a python script which took gene list outputs from a NGS experiment as input (either in ENSEMBL ID format or Conventional Gene name format) and returned the number of articles written on pubmed for each gene plus a user-inputted search term. For example, if the user input the search term "tendon development" or "ribosomal RNA", this would be searched with each gene in the user input gene list on the Pubmed API. Please see the EE282 writeup for more information.

For the for the updated program, I decided start by finding all orthologous genes in different species for every gene in the gene list, converting the ENSEMBL IDs into conventional gene names, and then performing the search against the Pubmed database via the [Pubmed API](https://www.ncbi.nlm.nih.gov/books/NBK25501/) for each gene and ortholog with the user inputted search terms.

The script has been converted into a command line program, called orthologue-finder-with-database.py, using the argparse python module, and takes 3 arguments as input:   

1) The file name of the excel file containing the input gene list in ENSEMBL ID format. The excel file should be formatted as shown in the following screenshot (Additional examples of input can be seen in the Sample Input Folder on this Github Repository):

![Sample Input](/Screenshots/Cartilage-Input.PNG)

2) The search term(s) you wish to search with each gene against the Pubmed API/Pubmed Research Database

3) The desired name of your output excel file, e.g. "output-file.xlsx"

Additionally, the program currently requires that you have downloaded the "ensembl_gene_lists" folder (and all files within), as well as the "orthologous_gene_lists" folder (and all files within). After entering the arguments for the program, it will ask for the paths to these folders, and the path to the parent folder containing your input file as well. Once these have all been input, the program should run normally, and output a file similar to the following screenshot (Additional examples of output can be seen in the Sample Output Folder on this Github Repository):

![Sample Output 1](/Screenshots/Cartilage-Output-1.PNG)

![Sample Output 2](/Screenshots/Cartilage-Output-2.PNG)

and we can validate that the number of articles written for this particular search (the search term in this case is "cartilage") by manually searching a gene and the search term on Pubmed, and looking at the number of articles returned on the web page.

Example 1, compared with genes boxed on first output screenshot:
![Validation 1](/Screenshots/Validation-1.PNG)

Example 2:
![Validation 2](/Screenshots/Validation-2.PNG)

The reason that these files are needed locally are because the program reads off of the correct orthologous gene list (downloaded off of [Biomart](https://www.ensembl.org/biomart/martview/a2d7e8d51d4d452d497f784ce4a5be89) to find all orthologous for each gene in the user-input gene list, and uses the correct ensembl gene list to translate ENSEMBL IDs into conventional gene names for proper searching against the Pubmed API. In the EE282 final project script, the ENSEMBL ID translation was done via the [ENSEMBL REST API](https://rest.ensembl.org/). This was originally considered for the first iteration for this project, but I found that when you added hundreds or thousands of genes (especially when adding in orthologous genes for multiple species), the API would restrict the number of API calls coming from the program (1x GET request per ENSEMBL ID), drastically slowing down the run time, or crashing the program entirely. The functions for the ENSEMBL API calls are still included within the program, although they are not being used currently, since local files are used instead. If interested, the functions for these calls are the "homolog_finder" and "multi_species_homfinder" functions for finding homologs/orthologous genes for one ENSEMBL ID for one species, and the "ens_id_finder" function for translating ENSEMBL IDs to conventional gene name.  The temporary solution was to create a local data storage/dictionary of sorts in the form of text files. A future improvement would be to only include API calls to the ENSEMBL REST API when the ENSEMBL ID was not found in the local ensembl ID dictionary/data storage text file.


The workflow for the program is as follows:

1) The user input excel file is read in as a pandas dataframe, and the species of the input ENSEMBL IDs is detected via the first 3-4 letters after the "ENS" prefix in the ENSEMBL ID. Currently, only zebrafish, human, rat, and mouse gene ENSEMBL IDs are supported (i.e. any of these species' gene IDs can be interconverted and the correct conventional gene name can be identified). Other Species may be added in future. The correct "local_multi_species_homfinder" function is used to find all orthologues for each gene in the input gene list excel file for each of the 3 other organisms (i.e. if the input gene list is human ENSEMBL IDs, then the function will find all orthologues for each ENSEMBL ID for zebrafish, rat, and mouse).

2) The dataframe is the restructured using the restruc_df() function to add 2 empty rows under each row, to allow for input of future information (the conventional gene name and the number of Pubmed articles)

3) The local_ens_id_finder() function looks at the correct ensembl gene list text file and converts every ENSEMBL ID (from the input gene or from any of the orthologues added in step 1) in the dataframe and adds them in the cell underneath their respective ENSEMBL ID

4) the pubcrawl function runs through each conventional gene name, and searches it along with the user inputted search term(s) against the Entrez eutils Pubmed API, and uses the bs4 package to parse out the number of articles written from the returned XML/HTML data, and place it in the cell under the corresponding conventional gene name in the dataframe.

5) A column is added at the end of the dataframe to calculate the total number of articles written for all orthologues for all species for each gene in the original input gene list. The dataframe is then sorted in descending order by the Totals column, and is then output as an excel file in the current directory according to the user-specified output excel file name.


Screenshots of Program Working:


After running program in command line, entering all arguments and paths correctly, the program should display a progress bar via the tqdm module, shown below:

![Working-Program](/Screenshots/program_working.PNG)

The Input file is in the correct path, shown below:

![program-input](/Screenshots/program_input.PNG)

After completion, the program should display a completed progress bar and the time taken to complete the program, shown below:

![Completed-Program](/Screenshots/program_completed.PNG)

The output excel file should now have been created in the working directory:

![program-output](/Screenshots/program_output.PNG)

## Known Issues and Potential Future Directions

1) When Genes in the user inputted excel file are not found in the local data storage text file, the output is returned under that gene as GENE NOT FOUND IN DATABASE, see example below:

![Not Found](/Screenshots/Genes_not_found.PNG)

An improvement could be made here where the program would use ENSEMBLE REST API requests only for genes not found in the local database, so that potentially useful genes are not missed in the search.

2) The program only takes ENSEMBL IDs as input currently, and only finds orthologous between zebrafish, human, rat and mouse. Interconversions between different gene ID formats can be considered in future, as well as the addition of orthologous genes from other organisms (e.g. Drosophila, C. Elegans, etc.)

3) Sometimes, the number of articles returned does not correspond to the number of articles written if you were to validate against the Pubmed website, and this is due to differences in what the Entrez Eutils API returns versus what the Pubmed website returns. However, I usually see that even if there is a discrepancy in the quantity of articles written, I have not seen a discrepancy where the API returns 0 articles, and the website returns a positive number of articles, so the program should still give some sort of useful information regardless of if the exact number of articles is 100% accurate.

4) Occasionally, if a search that returns 0 results on the Pubmed website will return an absurdly large number of articles (in the thousands to tens of thousands) in the program. This is again due to an issue with the Entrez Eutils API. A filter could be written to remove these abnormal results in future.

5) Despite being useful, the program could be more user friendly if turned into a web application with a User Interface to allow users who are not familiar with command line/programming usage to extract useful genes from their data.










