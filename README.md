# orthologue_finder_local_with_pubmed_search
Python program for automating literature searches for gene lists produced from Next Generation Sequencing (NGS) screens

import requests, sys
import pandas as pd
import bs4
import time
import numpy as np
from tqdm import tqdm



#from requests.adapters import HTTPAdapter
#from requests.packages.urllib3.util.retry import Retry
#def requests_retry_session(
#    retries=10,
#    backoff_factor=0.3,
#    status_forcelist=(500, 502, 504),
#    session=None,
#):
#    session = session or requests.Session()
#    retry = Retry(
#        total=retries,
#        read=retries,
#        connect=retries,
#        backoff_factor=backoff_factor,
#        status_forcelist=status_forcelist,
#    )
#    adapter = HTTPAdapter(max_retries=retry)
#    session.mount('http://', adapter)
#    session.mount('https://', adapter)
#    return session

startTime = time.time()
#IMPORTANT change this next line to the name of your ENSEMBL gene ID list. this is your input
orth_df = pd.read_excel('ens-gene-list-top-500.xlsx')

#This function finds all orthologs/homologs for an ensemble gene ID
# and puts them in a list together
# function inputs are ensembl gene id, and species for homolog needed
def homolog_finder(ens_gene_id, inputspecies):
    
    species = inputspecies


    server = "https://rest.ensembl.org"
    ext = "/homology/id/" + ens_gene_id + "?type=orthologues;format=condensed;target_species=" + species
     
    r = requests.get(server+ext, headers={ "Content-Type" : "text/xml", "User-Agent": "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.8) Gecko/20100722 Firefox/3.6.8 GTB7.1 (.NET CLR 3.5.30729)"})
 
     
 #need to create a way to re-try the http:get request if an error is returned.
 #if not r.ok:

    #sys.exit()
     
     
    #print(r.text)
    
    parse = bs4.BeautifulSoup(r.text)
    parse_hom = parse.select('homologies')


    homlist = []
    for x in range(0, len(parse_hom)):
        
        homlist.append( parse_hom[x]['id'])
    
    
    return homlist

#the below script reads an ensembl gene ID list and finds mouse homologs for all IDs, and places them in a table
#next to the original ID.


def  multi_species_homfinder(species_input):
    
    randgen = list(range(0, 101))
    counter = 1
    orth_df.insert(len(orth_df.columns), species_input + ' homologs', np.nan)
    homcount =  len(orth_df.columns)
    #iterate through gene list and run homolog_finder function to find gene homologs and place them in separate columns next to original Ensembl ID
    for i in tqdm(range(0, len(orth_df))):
        
        gentemp = homolog_finder(orth_df.iloc[i][0], species_input)
        
        #if there are multiple homologs to the gene for that species,
        #create new columns, and enter those into the new columns
        if len(gentemp) == 1:
            
            orth_df.iloc[i, homcount-1] = gentemp[0]
        
        
        #we want to check if the number of homologs is greater than the number of columns created for this species
        #if it is, then make enough new columns to fit in the total number of orthologues
        #If enough columns are made already, add the extra homologs into the existing columns    
          
        if len(gentemp) > counter:
                     
            for k in range(0, len(gentemp)-counter):
                
                orth_df.insert(homcount + k, species_input + ' Homologs' + ' ' + str(randgen[0]), np.nan)
                del randgen[0]
                counter = len(gentemp)
            for j in range(0, len(gentemp)):        
                
                orth_df.iloc[i, homcount-1 + j] = gentemp[j]
        
        #if no homologues exist, move to the next gene
        if len(gentemp) == 0:
            continue
        # if enough columns exist for the number of homologues, just add them into the existing columns sequentially
        elif len(gentemp) < counter and not len(gentemp) == 0:
           
            for m in range(0, len(gentemp)):
                
                orth_df.iloc[i, homcount-1 + m] = gentemp[m]
    
    
    return orth_df
        
#run the orthologue finder for mouse, zebrafish, human, and rat (recommended to choose only 1 or 2 homologs at most, otherwise the ensembl API may cut you off partway through because of too many HTTP requests)       
multi_species_homfinder('human')
#multi_species_homfinder('mouse')
#multi_species_homfinder('zebrafish')
#
#multi_species_homfinder('rat')
#multi_species_homfinder('drosophila')
#orth_df.to_excel('trial_homolog_finder.xlsx')

test_df = orth_df



#function to insert blank rows between all values of dataframe, to pre-allocate rows
#for finding gene names and pubmed counts
def pir(df):
    nans = np.where(np.empty_like(df.values), np.nan, np.nan)
    data = np.hstack([nans, df.values]).reshape(-1, df.shape[1])
    return pd.DataFrame(data, columns=df.columns)


#inserting blank rows, restructuring dataframe to fit future data
test_df = pir(pir(test_df))
test_df = test_df.append(test_df.iloc[0:3])
test_df = test_df.reset_index()
test_df = test_df.drop(test_df.columns[0], axis = 1)
test_df = test_df.drop([0,1,2], axis=0)
test_df = test_df.reset_index()
test_df = test_df.drop(test_df.columns[0], axis = 1)

#remove the extra NaNs so that there is only 2 NaN's after every orthologue
#after this, we can start filling in the NaN's with the actual gene names
for y in tqdm(range(3, len(test_df), 3)):
    try:
        test_df = test_df.drop(test_df.index[y])
    except IndexError:
        break

#reset the index    
test_df = test_df.reset_index()
test_df = test_df.drop(test_df.columns[0], axis = 1)


#this function finds the conventional gene name for all the ensembl IDs and
# places them directly 1 cell under the ensembl ID
def ens_id_finder(id_table):
    
    #start looping through the dataframe, by every 3, this is how the Ensembl IDs are organized
    for a in tqdm(range(0, len(id_table), 3)):
        for b in range(0, len(id_table.columns)):
            try:
                if type(id_table.iloc[a][b]) == str:
                #access ensembl API and look up each ensembl ID on the server
                    server = "http://rest.ensembl.org"
                    ext = "/lookup/id/" + id_table.iloc[a][b]

                #get the gene info for the ensembl ID
                    r = requests.get(server+ext, headers={ "Content-Type" : "application/json", "User-Agent": "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.8) Gecko/20100722 Firefox/3.6.8 GTB7.1 (.NET CLR 3.5.30729)"})
                         
                    #    if not r.ok:
                    #      r.raise_for_status()
                    #      sys.exit()
                    #

                #convert to json file
                    decoded = r.json()
                #if the ensembl ID is not found or is outdated, print error under the respective ensembl ID
                    if 'error' in decoded:
                        ens_gene_name = 'ERROR:NO ENSEMBL ID FOUND'
                        id_table.iloc[a+1][b] = ens_gene_name
                #else, put the correct gene name in the respective spot in the conventional gene list
                    else:
                        ens_gene_name = decoded['display_name']
                        id_table.iloc[a+1][b] = ens_gene_name
                else:
                    continue
            except IndexError:
                break
                
        
    return id_table

ens_id_finder(test_df)

#this function finds the number of pubmed articles in for each gene in the dataframe + your search term keyword

def pubcrawl(dataframe, keyword):
    #loop through the dataframe every 3 rows, and for every column
    for x in tqdm(range(1, len(dataframe), 3)):
        for y in range(0 ,len(dataframe.columns)):
            #if the gene name is NaN or a 0, then continue on
            if pd.isnull(dataframe.iloc[x,y]) or dataframe.iloc[x,y] == 0:
                continue
            #search the entrez eutility pubmed API for each gene plus your keyword
            res = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=' + dataframe.iloc[x, y] + '+'+ keyword)
            ncbi = bs4.BeautifulSoup(res.text)
            count = ncbi.select('Count')
            
            #if nothing is found or there is an error, simply make the number of articles as 0, and continue
            if ncbi.select('errorlist') != []:
                dataframe.iloc[x+1, y] = 0
                continue
            #if nothing is found, make the number of articles as 0 and continue
            if count == []:
                dataframe.iloc[x+1,y] = 0
                continue
            #once you have the number of papers for your search term, conver to string, remove the html tags, and convert to integer
            #then add it to the cell directly underneath the respective gene in the dataframe
            once = str(count[0])
            once = once.replace('<count>', '')
            once = once.replace('</count>', '')
            dataframe.iloc[x+1,y] = int(once)
    
#IMPORTANT, CHANGE THE SEARCH TERM, WRITTEN HERE AS TENDON, TO WHATEVER YOU WISH TO SEARCH WITH YOUR GENES ON PUBMED    
pubcrawl(test_df, 'mechanical')

#create a last column for totals
test_df = test_df.assign(Total= np.nan)

#sum up counts for all genes at the end and place in the last column for each row
for p in range(2, len(test_df), 3):
    test_df.iloc[p, len(test_df.columns)-1] = test_df.iloc[p, :].sum(axis=0)
    

#add counts to the end of every row to help sort
test_df.Total = test_df.Total.bfill()
#create a helper column and fill with numbers to help sort but keep dataframe formatting
test_df["helper"] = np.arange(len(test_df))//3
#sort by helper values in descending order (highest to lowest pubmed article counts)
test_df = test_df.sort_values(['Total', 'helper'], ascending = False)
#remove the helper column
test_df = test_df.drop(columns = "helper")


#remove excess counts (not necessary, can improve aesthetically if desired)
#for p in range(2, len(test_df), 3):
#    test_df.iloc[p-1, len(test_df.columns)-1] = np.nan
#    test_df.iloc[p-2, len(test_df.columns)-1] = np.nan
    
print ('The script took {0} seconds !'.format(time.time() - startTime))

#the output dataframe is called test_df
#if you want, you can output the dataframe to an excel file using the next line of code
test_df.to_excel('ens-id-output-top-500-mechanical.xlsx')


