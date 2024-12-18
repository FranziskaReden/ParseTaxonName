import pandas as pd
import sys
import os.path
import shutil
import string
import urllib
import urllib.request
import pathlib
from datetime import datetime

import utils

def sort_taxa_names(folder:pathlib.Path) -> pd.DataFrame: 

    file_name = os.path.join(os.path.join(folder, 'taxdmp'), 'names.dmp')

    print('Reading in file '+str(file_name)+'...')
    # Read in file line after line
    with open(file_name) as t: 
        taxa = t.readlines()

    # Split lines into corresponding columns
    for i in range (len(taxa)): 
        taxa[i] = utils.read_line(taxa[i])

    # Declare header for DataFrame
    header = ['tax_id', 'name_txt', 'unique name', 'name class']
    
    # Sort list of lists according to the second position in list (taxon name)
    print('Sorting taxa names alphabetically...')
    taxa.sort(key = lambda row: row[1].lower())

    # Transform sorted list of lists into DataFrame and write 
    taxa_df = pd.DataFrame(taxa, columns=header)
    taxa_df.to_csv(os.path.join(folder, 'taxa_names_sorted.tsv'), sep='\t', index=False)
    print(f'The DataFrame containing the sorted taxon names and taxon IDs were written into file: {os.path.join(folder, 'taxa_names_sorted.tsv')}.\n')

    return taxa_df

def get_indeces(folder:pathlib.Path, taxa:pd.DataFrame) -> dict:
    
    print('Creating lexicon for taxa names DataFrame...')
    # Create dictionary that will store the starting and ending indices for each letter in the alphabet
    indeces = {}

    # Get list of letters in alphabet
    letters = string.ascii_uppercase
    for letter in letters: 
        indeces[letter] = []    

    # Get indices in the taxa df for all names not starting with a letter. Starting index is 0.
    indeces.setdefault('_', []).append(0)
    for i in range(len(taxa.index)): 
        if taxa['name_txt'][i][0].upper() == letters[0]: 
            indeces.setdefault('_', []).append(i-1)
            break
    
    # For each letter get the starting and ending index (i) in the taxa DataFrame
    for key in letters: 
        indeces.setdefault(key, []).append(i)
        while i < len(taxa.index): 
            name = taxa['name_txt'][i][0].upper()
            # Once the first letter in the name does not correspond to the current letter anymore, store index in indeces dict and break
            if name != key:
                indeces.setdefault(key, []).append(i-1)
                break 
            i += 1
        if i == (len(taxa.index)): 
            indeces.setdefault(key, []).append(i-1)

    # Write results into taxa_indeces.txt file
    with open(os.path.join(folder, 'taxa_indeces.txt'), 'w') as w: 
        for key in indeces: 
            w.write(str(key)+' '+str(indeces[key][0])+' '+str(indeces[key][1])+'\n')
        print('Lexicon for taxon name DataFrame was written into file: '+str(os.path.join(folder, 'taxa_indeces.txt'))+'\n')

    return indeces

def read_indices(file:pathlib.Path) -> dict: 
    
    with open (file) as t: 
        list_index_raw = t.readlines()
        # Write indeces into a dictionary
        list_index = {}
        for i in range (len(list_index_raw)): 
            list_index_raw[i] = list_index_raw[i].strip('\n')
            tmp = list_index_raw[i].split(' ')
            list_index[tmp[0]] = [int(tmp[1]), int(tmp[2])]

    return list_index

def get_nodes_file(folder:pathlib.Path, names:pd.DataFrame) -> pd.DataFrame: 

    nodes_file = os.path.join(os.path.join(folder, 'taxdmp'), 'nodes.dmp')

    print('Reading in file '+str(nodes_file)+'...')
    # Read in file line after line
    with open(nodes_file) as t: 
        nodes = t.readlines()

    for i in range (len(nodes)): 
        nodes[i] = utils.read_line(nodes[i])        

    # Declare header for DataFrame
    header = ['tax_id', 'parent_tax_id', 'rank', 'name class', 'embl code', 'division id', \
              'inherited div flag', 'genetic code id', 'inherited GC  flag', 'mitochondrial genetic code id', \
                'inherited MGC flag', 'GenBank hidden flag', 'hidden subtree root flag'] 

    taxa_df = names[['tax_id', 'name_txt']]

    nodes_df = pd.DataFrame(nodes, columns=header)
    nodes_df = nodes_df[['tax_id', 'parent_tax_id', 'rank']]
    nodes_df = nodes_df.astype({'tax_id': int, 'parent_tax_id': int})
    taxa_df = taxa_df.astype({'tax_id': int})

    nodes_df = nodes_df.merge(taxa_df, how='left', on='tax_id')
    nodes_df.to_csv(os.path.join(folder, 'nodes.tsv'), sep='\t', index=False)
    print('The DataFrame containing the the lineages was written into: '+os.path.join(folder, 'nodes.tsv')+'\n')

    nodes_df = pd.read_csv(os.path.join(folder, 'nodes.tsv'), sep='\t')
    nodes_df.head()

    return nodes_df

def get_dumpfile(folder, timeout=540): 

    print('Downloading taxdmp.zip file from the NCBI taxonomy database...\n')
    url = 'https://ftp.ncbi.nih.gov/pub/taxonomy/taxdmp.zip'

    # Download names.dmp file
    try: 
        with urllib.request.urlopen(url, timeout=timeout) as dl_file:
            with open(os.path.join(folder,'taxdmp.zip'), 'wb') as out_file:
                out_file.write(dl_file.read())

    except Exception as x: 
        print(x)
        sys.exit(1)

    date = datetime.now()
    with open(os.path.join(folder, 'update.log'), 'a') as w: 
        w.write(f'NCBI taxonomy last downloaded and updated on: {date}.')

def get_taxa(folder:str) -> list: 

    if os.path.isdir(folder) is False: 
        os.mkdir(folder)

    if os.path.exists(os.path.join(folder, 'taxdmp.zip')) is False: 
        get_dumpfile(folder)  

    if os.path.exists(os.path.join(folder, 'taxdmp')) is False: 
        shutil.unpack_archive(filename=os.path.join(folder, 'taxdmp.zip'), extract_dir=os.path.join(folder, 'taxdmp'))

    if os.path.exists(os.path.join(folder,'taxa_names_sorted.tsv')) is False: 
        # Sort taxa names; get indices to create lexicon.
        taxa = sort_taxa_names(folder)
        list_index = get_indeces(folder, taxa)
    else: 
        print('Reading in '+os.path.join(folder, 'taxa_names_sorted.tsv')+' file...')
        taxa = pd.read_csv(os.path.join(folder, 'taxa_names_sorted.tsv'), sep='\t')
        list_index = read_indices(os.path.join(folder, 'taxa_indeces.txt'))

    return taxa, list_index

def get_nodes(folder_name:str) -> pd.DataFrame: 

    folder = folder_name

    if os.path.exists(os.path.join(folder, 'nodes.tsv')) is False: 
        taxa, list_index = get_taxa(folder)
        nodes_df = get_nodes_file(folder, taxa[taxa['name class'] == 'scientific name'])
    else: 
        print('Reading in '+os.path.join(folder, 'nodes.tsv')+' file...')
        nodes_df = pd.read_csv(os.path.join(folder, 'nodes.tsv'), sep='\t')

    nodes_df = nodes_df.set_index('tax_id')

    return nodes_df

def update_db(folder): 

    if os.path.isdir(folder) is False: 
        os.mkdir(folder)

    get_dumpfile(folder) 
    shutil.unpack_archive(filename=os.path.join(folder, 'taxdmp.zip'), extract_dir=os.path.join(folder, 'taxdmp'))
    taxa = sort_taxa_names(folder)
    list_index = get_indeces(folder, taxa)
    nodes_df = get_nodes_file(folder, taxa[taxa['name class'] == 'scientific name'])