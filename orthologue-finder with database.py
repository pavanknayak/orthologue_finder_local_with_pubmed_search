import requests
# import sys
import pandas as pd
import bs4
import time
import numpy as np
from tqdm import tqdm
import os

gene_list_path = r"C:\Users\Pavan Nayak\Documents\Programming Practice\Tools\ensembl_gene_lists"
orth_list_path = r"C:\Users\Pavan Nayak\Documents\Programming Practice\Tools\orthologous_gene_lists"
input_file_path = r"C:\Users\Pavan Nayak\Documents\Programming Practice\Tools"
human_gl = os.path.join(gene_list_path, "human_genes.txt")
zebrafish_gl = os.path.join(gene_list_path, "zebrafish_genes.txt")
mouse_gl = os.path.join(gene_list_path, "mouse_genes.txt")
rat_gl = os.path.join(gene_list_path, "rat_genes.txt")

zeb2hum = os.path.join(orth_list_path, "zebrafish2human.txt")
zeb2mou = os.path.join(orth_list_path, "zebrafish2mouse.txt")
zeb2rat = os.path.join(orth_list_path, "zebrafish2rat.txt")
hum2zeb = os.path.join(orth_list_path, "human2zebrafish.txt")
hum2mou = os.path.join(orth_list_path, "human2mouse.txt")
hum2rat = os.path.join(orth_list_path, "human2rat.txt")
mou2hum = os.path.join(orth_list_path, "mouse2human.txt")
mou2zeb = os.path.join(orth_list_path, "mouse2zebrafish.txt")
mou2rat = os.path.join(orth_list_path, "mouse2rat.txt")
rat2hum = os.path.join(orth_list_path, "rat2human.txt")
rat2zeb = os.path.join(orth_list_path, "rat2zebrafish.txt")
rat2mou = os.path.join(orth_list_path, "rat2mouse.txt")

human_genelist = pd.read_csv(human_gl)
zebrafish_genelist = pd.read_csv(zebrafish_gl)
mouse_genelist = pd.read_csv(mouse_gl)
rat_genelist = pd.read_csv(rat_gl)

zebrafish2human = pd.read_csv(zeb2hum)
zebrafish2mouse = pd.read_csv(zeb2mou)
zebrafish2rat = pd.read_csv(zeb2rat)
human2zebrafish = pd.read_csv(hum2zeb)
human2mouse = pd.read_csv(hum2mou)
human2rat = pd.read_csv(hum2rat)
mouse2human = pd.read_csv(mou2hum)
mouse2rat = pd.read_csv(mou2rat)
mouse2zebrafish = pd.read_csv(mou2zeb)
rat2human = pd.read_csv(rat2hum)
rat2zebrafish = pd.read_csv(rat2zeb)
rat2mouse = pd.read_csv(rat2mou)

startTime = time.time()
# IMPORTANT change this next line to the name of your ENSEMBL gene ID list. this is your input
orth_df = pd.read_excel(os.path.join(input_file_path, 'test-ens-gene-list-short1.xlsx'))


# these four functions are to search the dataframes for homologous genes for a particular gene input for each different species
# I created a local gene ortho/homology finder for each species for human, mouse, rat, and zebrafish so far
# each function takes as input a string ensembl ID and returns 3 lists, one list for the homologs for each species for your input gene
# These will be called in the local multispecies homolog finder created later to change the overall input dataframe into one containing all the homologs for all genes in the list

def local_homfinder_human(ens_string):
    hl_human2rat = []
    hl_human2mouse = []
    hl_human2zebrafish = []

    hom_indexlist1 = list(np.where(human2rat["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_human2rat.append(human2rat.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(human2mouse["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_human2mouse.append(human2mouse.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(human2zebrafish["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_human2zebrafish.append(human2zebrafish.iloc[hom_indexlist3[c], 1])

    return hl_human2rat, hl_human2mouse, hl_human2zebrafish


def local_homfinder_rat(ens_string):
    hl_rat2human = []
    hl_rat2mouse = []
    hl_rat2zebrafish = []

    hom_indexlist1 = list(np.where(rat2human["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_rat2human.append(rat2human.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(rat2mouse["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_rat2mouse.append(rat2mouse.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(rat2zebrafish["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_rat2zebrafish.append(rat2zebrafish.iloc[hom_indexlist3[c], 1])

    return hl_rat2human, hl_rat2mouse, hl_rat2zebrafish


def local_homfinder_mouse(ens_string):
    hl_mouse2human = []
    hl_mouse2rat = []
    hl_mouse2zebrafish = []

    hom_indexlist1 = list(np.where(mouse2human["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_mouse2human.append(mouse2human.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(mouse2rat["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_mouse2rat.append(mouse2rat.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(mouse2zebrafish["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_mouse2zebrafish.append(mouse2zebrafish.iloc[hom_indexlist3[c], 1])

    return hl_mouse2human, hl_mouse2rat, hl_mouse2zebrafish


def local_homfinder_zebrafish(ens_string):
    hl_zebrafish2human = []
    hl_zebrafish2rat = []
    hl_zebrafish2mouse = []

    hom_indexlist1 = list(np.where(zebrafish2human["Gene stable ID"] == ens_string)[0])
    for a in range(0, len(hom_indexlist1)):
        hl_zebrafish2human.append(zebrafish2human.iloc[hom_indexlist1[a], 1])

    hom_indexlist2 = list(np.where(zebrafish2rat["Gene stable ID"] == ens_string)[0])
    for b in range(0, len(hom_indexlist2)):
        hl_zebrafish2rat.append(zebrafish2rat.iloc[hom_indexlist2[b], 1])

    hom_indexlist3 = list(np.where(zebrafish2mouse["Gene stable ID"] == ens_string)[0])
    for c in range(0, len(hom_indexlist3)):
        hl_zebrafish2mouse.append(zebrafish2mouse.iloc[hom_indexlist3[c], 1])

    return hl_zebrafish2human, hl_zebrafish2rat, hl_zebrafish2mouse

#These four functions call the previous local_homfinder functions and use them to identify all of the homologs for all genes in an
#inputted dataframe.
def local_multi_species_homfinder_zeb(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["human", "rat", "mouse"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_zebrafish(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df

def local_multi_species_homfinder_hum(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["rat", "mouse", "zebrafish"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_human(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df


def local_multi_species_homfinder_rat(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["human", "mouse", "zebrafish"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_rat(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df


def local_multi_species_homfinder_mouse(orth_df):
    for a in range(0, 3):

        counter = 1
        namegen = list(range(0, 500))
        namelist = ["human", "rat", "zebrafish"]
        homcount = len(orth_df.columns)
        for b in range(0, len(orth_df)):
            namecount = len(orth_df.columns)
            genetemp = local_homfinder_mouse(orth_df.iloc[b, 0])
            if len(genetemp[a]) == 1 and namecount == homcount:
                orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) == 1 and not namecount == homcount:
                orth_df.iloc[b, homcount] = genetemp[a][0]

            elif len(genetemp[a]) > counter:

                for c in range(0, len(genetemp[a]) - counter + 1):
                    orth_df.insert(namecount + c, namelist[a] + ' homolog ' + str(namegen[0]), np.nan)
                    del namegen[0]
                    counter = len(genetemp[a])
                for d in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + d] = genetemp[a][d]

            if len(genetemp[a]) == 0:
                continue

            elif len(genetemp[a]) < counter and not len(genetemp[a]) == 0:

                for e in range(0, len(genetemp[a])):
                    orth_df.iloc[b, homcount + e] = genetemp[a][e]

    return orth_df


#function to insert blank rows between all values of dataframe, to pre-allocate rows
#for finding gene names and pubmed counts
def pir(df):
    nans = np.where(np.empty_like(df.values), np.nan, np.nan)
    data = np.hstack([nans, df.values]).reshape(-1, df.shape[1])
    return pd.DataFrame(data, columns=df.columns)


#function for inserting blank rows, restructuring dataframe to fit future data
def restruc_df(test_df):
    test_df = pir(pir(test_df))
    test_df = test_df.append(test_df.iloc[0:3])
    test_df = test_df.reset_index()
    test_df = test_df.drop(test_df.columns[0], axis = 1)
    test_df = test_df.drop([0,1,2], axis=0)
    test_df = test_df.reset_index()
    test_df = test_df.drop(test_df.columns[0], axis = 1)

    #remove the extra NaNs so that there is only 2 NaN's after every orthologue
    #after this, we can start filling in the NaN's with the actual gene names
    for y in range(3, len(test_df), 3):
        try:
            test_df = test_df.drop(test_df.index[y])
        except IndexError:
            break

    #reset the index
    test_df = test_df.reset_index()
    test_df = test_df.drop(test_df.columns[0], axis = 1)

    return test_df

#local ensembl ID finder. After restructuring the dataframe to add spaces between each row item, we call this function
#to do find the actual gene name for each ensembl ID in the dataframe, and add it to the cell directly underneath the respective ensembl_ID
#then output the dataframe.
def local_ens_id_finder(id_table):
    for a in range(0, len(id_table), 3):
        for b in range(0, len(id_table.columns)):
            if type(id_table.iloc[a][b]) == str:

                if id_table.iloc[a][b].startswith('ENSDARG'):

                    genindex = list(np.where(zebrafish_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = zebrafish_genelist.iloc[genindex[0]][1]


                elif id_table.iloc[a][b].startswith('ENSG'):

                    genindex = list(np.where(human_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = human_genelist.iloc[genindex[0]][1]

                elif id_table.iloc[a][b].startswith('ENSRNOG'):

                    genindex = list(np.where(rat_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = rat_genelist.iloc[genindex[0]][1]

                elif id_table.iloc[a][b].startswith('ENSMUSG'):

                    genindex = list(np.where(mouse_genelist["Gene stable ID"] == id_table.iloc[a][b])[0])
                    if len(genindex) == 0:
                        id_table.iloc[a + 1][b] = "GENE NOT FOUND IN DATABASE"
                    else:
                        id_table.iloc[a + 1][b] = mouse_genelist.iloc[genindex[0]][1]

    return id_table
#EVERYTHING FROM THIS POINT ON WAS THE ORIGINAL SCRIPT CALLING APIs. Will integrate into overall program soon
#This function finds all orthologs/homologs for an ensemble gene ID
# and puts them in a list together USING THE ENSEMBL API
# function inputs are ensembl gene id, and species for homolog needed
def homolog_finder(ens_gene_id, inputspecies):
    species = inputspecies
    server = "https://rest.ensembl.org"
    ext = "/homology/id/" + ens_gene_id + "?type=orthologues;format=condensed;target_species=" + species
    r = requests.get(server+ext, headers={ "Content-Type" : "text/xml", "User-Agent": "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.8) Gecko/20100722 Firefox/3.6.8 GTB7.1 (.NET CLR 3.5.30729)"})
    parse = bs4.BeautifulSoup(r.text)
    parse_hom = parse.select('homologies')


    homlist = []
    for x in range(0, len(parse_hom)):
        
        homlist.append( parse_hom[x]['id'])
    return homlist

#the below function reads an ensembl gene ID list and finds mouse homologs for all IDs, and places them in a table
#next to the original ID, again USING THE ENSEMBL ID


def  multi_species_homfinder(species_input):
    
    randgen = list(range(0, 101))
    counter = 1
    orth_df.insert(len(orth_df.columns), species_input + ' homologs', np.nan)
    homcount =  len(orth_df.columns)
    #iterate through gene list and run homolog_finder function to find gene homologs and place them in separate columns next to original Ensembl ID
    for i in tqdm(range(0, len(orth_df))):
        namecount = len(orth_df.columns)
        gentemp = homolog_finder(orth_df.iloc[i][0], species_input)
        
        #if there are multiple homologs to the gene for that species,
        #create new columns, and enter those into the new columns
        if len(gentemp) == 1 and namecount == homcount:
            orth_df.insert(namecount, namelist[a] + ' homolog', np.nan)
            orth_df.iloc[i, homcount-1] = gentemp[0]

        elif len(genetemp[a]) == 1 and not namecount == homcount:
            orth_df.iloc[b, homcount] = genetemp[a][0]
        #we want to check if the number of homologs is greater than the number of columns created for this species
        #if it is, then make enough new columns to fit in the total number of orthologues
        #If enough columns are made already, add the extra homologs into the existing columns    
          
        elif len(gentemp) > counter:
                     
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
#multi_species_homfinder('human')
#multi_species_homfinder('mouse')
#multi_species_homfinder('zebrafish')
#multi_species_homfinder('rat')
#multi_species_homfinder('drosophila')
#orth_df.to_excel('trial_homolog_finder.xlsx')

#test_df = orth_df





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
    # create a last column for totals
    dataframe = dataframe.assign(Total= np.nan)
    # sum up counts for all genes at the end and place in the last column for each row
    for p in range(2, len(dataframe), 3):
        dataframe.iloc[p, len(dataframe.columns) - 1] = dataframe.iloc[p, :].sum(axis=0)

    # add counts to the end of every row to help sort
    dataframe.Total = dataframe.Total.bfill()
    # create a helper column and fill with numbers to help sort but keep dataframe formatting
    dataframe["helper"] = np.arange(len(dataframe))//3
    # sort by helper values in descending order (highest to lowest pubmed article counts)
    dataframe = dataframe.sort_values(['Total', 'helper'], ascending = False)
    # remove the helper column
    dataframe = dataframe.drop(columns = "helper")
    # remove excess counts (not necessary, can improve aesthetically if desired)
    for p in range(2, len(dataframe), 3):
        dataframe.iloc[p - 1, len(dataframe.columns) - 1] = np.nan
        dataframe.iloc[p - 2, len(dataframe.columns) - 1] = np.nan
    return dataframe
#IMPORTANT, CHANGE THE SEARCH TERM, WRITTEN HERE AS TENDON, TO WHATEVER YOU WISH TO SEARCH WITH YOUR GENES ON PUBMED    
#pubcrawl(test_df, 'mechanical')





    

#add counts to the end of every row to help sort
#test_df.Total = test_df.Total.bfill()
#create a helper column and fill with numbers to help sort but keep dataframe formatting
#test_df["helper"] = np.arange(len(test_df))//3
#sort by helper values in descending order (highest to lowest pubmed article counts)
#test_df = test_df.sort_values(['Total', 'helper'], ascending = False)
#remove the helper column
#test_df = test_df.drop(columns = "helper")



    
print ('The script took {0} seconds !'.format(time.time() - startTime))

#the output dataframe is called test_df
#if you want, you can output the dataframe to an excel file using the next line of code
#test_df.to_excel('ens-id-output-top-500-mechanical.xlsx')
