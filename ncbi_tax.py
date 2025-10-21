import sys
import os.path
import shutil
import string
import urllib
import urllib.request
import pathlib
from datetime import datetime
import pandas as pd
from collections import defaultdict
import json
from tqdm import tqdm

import utils

def sort_taxa_names(folder:pathlib.Path) -> pd.DataFrame:

    file_name = os.path.join(os.path.join(folder, 'taxdmp'), 'names.dmp')

    print('Reading in file '+str(file_name)+'...')
    # Read in file line after line
    with open(file_name, encoding='utf=8') as t:
        taxa = t.readlines()

    # Split lines into corresponding columns
    for i, value in enumerate(taxa):
        taxa[i] = utils.read_line(taxa[i])

    # Declare header for DataFrame
    header = ['tax_id', 'name_txt', 'unique name', 'name class']

    # Sort list of lists according to the second position in list (taxon name)
    print('Sorting taxa names alphabetically...')
    taxa.sort(key = lambda row: row[1].lower())

    # Transform sorted list of lists into DataFrame and write
    taxa_df = pd.DataFrame(taxa, columns=header)
    taxa_df.to_csv(os.path.join(folder, 'taxa_names_sorted.tsv'), sep='\t', index=False)
    print('The DataFrame containing the sorted taxon names and taxon IDs \
were written into file: '+os.path.join(folder, 'taxa_names_sorted.tsv')+'.\n')

    return taxa_df

def get_indeces(folder:pathlib.Path, taxa:pd.DataFrame) -> dict:

    print('Creating lexicon for taxa names DataFrame...')
    # Create dictionary that will store the starting and
    # ending indices for each letter in the alphabet
    indeces = {}

    # Get list of letters in alphabet
    letters = string.ascii_uppercase
    for letter in letters:
        indeces[letter] = []

    # Get indices in the taxa df for all names not starting with a letter.
    # Starting index is 0.
    indeces.setdefault('_', []).append(0)
    for i in range(len(taxa.index)):
        if taxa.at[i, 'name_txt'][0].upper() == letters[0]:
            indeces.setdefault('_', []).append(i-1)
            break

    # For each letter get the starting and ending index (i) in the taxa DataFrame
    for key in letters:
        indeces.setdefault(key, []).append(i)
        while i < len(taxa.index):
            name = taxa.at[i,'name_txt'][0].upper()
            # Once the first letter in the name does not correspond
            # to the current letter anymore, store index in indeces dict and break
            if name != key:
                indeces.setdefault(key, []).append(i-1)
                break
            i += 1
        if i == (len(taxa.index)):
            indeces.setdefault(key, []).append(i-1)

    # Write results into taxa_indeces.txt file
    with open(os.path.join(folder, 'taxa_indeces.txt'), 'w', encoding='utf-8') as w:
        for key, item in indeces.items():
            w.write(str(key)+' '+str(item[0])+' '+str(item[1])+'\n')
        print('Lexicon for taxon name DataFrame was written into file: \
              '+str(os.path.join(folder, 'taxa_indeces.txt'))+'\n')

    return indeces

def read_indices(file:pathlib.Path) -> dict:

    with open (file, encoding='utf=8') as t:
        list_index_raw = t.readlines()
        # Write indeces into a dictionary
        list_index = {}
        for i, value in enumerate(list_index_raw):
            list_index_raw[i] = list_index_raw[i].strip('\n')
            tmp = list_index_raw[i].split(' ')
            list_index[tmp[0]] = [int(tmp[1]), int(tmp[2])]

    return list_index

def get_nodes_file(folder:pathlib.Path, names:pd.DataFrame) -> pd.DataFrame:
    '''
    Function to get the nodes DataFrame from the NCBI taxonomy database.
    Writes nodes.tsv file.

    Parameters
    ----------
    folder : str
        Path to folder in which to write the nodes.tsv file.
    names : pd.DataFrame
        DataFrame holding the taxa names and tax IDs.

    Returns
    ----------
    nodes_df : pd.DataFrame
        DataFrame holding the nodes of the NCBI taxonomy database.
    '''

    nodes_file = os.path.join(os.path.join(folder, 'taxdmp'), 'nodes.dmp')

    print('Reading in file '+str(nodes_file)+'...')
    # Read in file line after line
    with open(nodes_file, encoding='utf=8') as t:
        nodes = t.readlines()

    for i, value in enumerate(nodes):
        nodes[i] = utils.read_line(nodes[i])

    # Declare header for DataFrame
    header = ['tax_id', 'parent_tax_id', 'rank', 'name class', 'embl code',
        'division id', 'inherited div flag', 'genetic code id',
        'inherited GC  flag', 'mitochondrial genetic code id', 
        'inherited MGC flag', 'GenBank hidden flag', 'hidden subtree root flag'] 

    taxa_df = names[['tax_id', 'name_txt']]

    nodes_df = pd.DataFrame(nodes, columns=header)
    nodes_df = nodes_df[['tax_id', 'parent_tax_id', 'rank']]
    nodes_df = nodes_df.astype({'tax_id': int, 'parent_tax_id': int})
    taxa_df = taxa_df.astype({'tax_id': int})

    nodes_df = nodes_df.merge(taxa_df, how='left', on='tax_id')
    nodes_df.to_csv(os.path.join(folder, 'nodes.tsv'), sep='\t', index=False)
    print('The DataFrame containing the the lineages was written into: \
          '+os.path.join(folder, 'nodes.tsv')+'\n')

    nodes_df = pd.read_csv(os.path.join(folder, 'nodes.tsv'), sep='\t')
    nodes_df.head()

    return nodes_df

def get_dumpfile(folder, timeout=540) -> None:
    '''
    Function to download the taxdmp.zip file from the NCBI taxonomy database.

    Parameters
    ----------
    folder : str
        Path to folder in which to write the taxdmp.zip file.
    timeout : int
        Timeout in seconds for downloading the file. Default is 540 seconds.
    '''

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
    with open(os.path.join(folder, 'update.log'), 'a', encoding='utf-8') as w:
        w.write(f'NCBI taxonomy last downloaded and updated on: {date}.\n')

def get_taxa(folder:str) -> list:
    '''
    Function to get the taxa DataFrame and the lexicon from the NCBI taxonomy database
    If taxa_names_sorted.tsv and taxa_indeces.txt files do not exist, they will be created.

    Parameters
    ----------
    folder : str
        Path to folder in which to find/write the NCBI taxonomy database files.

    Returns
    ----------
    taxa : pd.DataFrame
        DataFrame holding the taxa names and tax IDs.
    list_index : dict
        Dictionary holding the lexicon for the taxa DataFrame.
    '''

    if os.path.isdir(folder) is False:
        os.mkdir(folder)

    if os.path.exists(os.path.join(folder, 'taxdmp.zip')) is False:
        get_dumpfile(folder)

    if os.path.exists(os.path.join(folder, 'taxdmp')) is False:
        shutil.unpack_archive(filename=os.path.join(folder, 'taxdmp.zip'),
                              extract_dir=os.path.join(folder, 'taxdmp'))

    if os.path.exists(os.path.join(folder,'taxa_names_sorted.tsv')) is False:
        # Sort taxa names; get indices to create lexicon.
        taxa = sort_taxa_names(folder)
        list_index = get_indeces(folder, taxa)
    else:
        print('Reading in '+os.path.join(folder, 'taxa_names_sorted.tsv')+' file...')
        taxa = pd.read_csv(os.path.join(folder, 'taxa_names_sorted.tsv'), sep='\t')
        list_index = read_indices(os.path.join(folder, 'taxa_indeces.txt'))

    return taxa, list_index

def get_nodes(folder:str) -> pd.DataFrame:
    '''
    Function to get the nodes DataFrame from the NCBI taxonomy database.
    If nodes.tsv file does not exist, it will be created.

    Parameters
    ----------
    folder : str
        Path to folder in which to find/write the nodes.tsv file.

    Returns
    ----------
    nodes_df : pd.DataFrame
        DataFrame holding the nodes of the NCBI taxonomy database.
    '''

    if os.path.exists(os.path.join(folder, 'nodes.tsv')) is False:
        taxa, list_index = get_taxa(folder)
        nodes_df = get_nodes_file(folder, taxa[taxa['name class'] == 'scientific name'])
    else:
        print('Reading in '+os.path.join(folder, 'nodes.tsv')+' file...')
        nodes_df = pd.read_csv(os.path.join(folder, 'nodes.tsv'), sep='\t')

    nodes_df = nodes_df.set_index('tax_id')

    return nodes_df

def get_homonyms_file(folder:str, taxa:pd.DataFrame, redo = False) -> list:
    '''
    Function to get all homonyms in the taxa DataFrame.	
    Writes homonyms into a json file.

    Parameters
    ----------
    folder : str
        Path to folder in which to write the homonyms file.
    taxa : pd.DataFrame
        DataFrame holding the taxa names and tax IDs.
    '''

    homonyms_file = os.path.join(folder, 'homonyms.json')

    if os.path.exists(homonyms_file) and redo is False:
        with open(homonyms_file, "r") as file:
            duplicates_dict = json.load(file)
        return duplicates_dict

    homonyms = defaultdict(list)
    for idx, row in tqdm(taxa.iterrows(), total=len(taxa)): 
        homonyms[row['name_txt']].append(idx)
        # If we encounter a homonym, set dup column ccordingly
        if len(homonyms[row['name_txt']]) > 1:
            for idx in homonyms[row['name_txt']]:
                taxa.at[idx, 'dup'] = 1

    # keep only duplicates
    duplicates_dict = {name: ids for name, ids in homonyms.items() if len(ids) > 1}

    with open(homonyms_file, "w") as f:
        json.dump(duplicates_dict, f, indent=4)
    print('Homonyms were written into file '+homonyms_file+'.')

    return duplicates_dict

def add_dup_to_taxa(folder:str, taxa_df:pd.DataFrame):
    '''
    Function to add a 'dup' column to the taxa DataFrame.

    Parameters
    ----------
    folder : str
        Path to folder in which to write the updated taxa DataFrame.
    taxa_df : pd.DataFrame
        DataFrame holding the taxa names and tax IDs.
    '''

    if 'dup' in taxa_df.columns:
        return
    
    taxa_df['dup'] = 0
    _ = get_homonyms_file(folder, taxa_df, redo = True)       
    
    taxa_df.to_csv(os.path.join(folder, 'taxa_names_sorted.tsv'), sep='\t', index=False)
    return

def update_db(folder: str):
    '''
    Function to update the NCBI taxonomy database.

    Parameters
    ----------
    folder : str
        Path to folder in which to write/find the NCBI taxonomy database files.
    '''

    if os.path.isdir(folder) is False:
        os.mkdir(folder)

    get_dumpfile(folder)
    shutil.unpack_archive(filename=os.path.join(folder, 'taxdmp.zip'),
                          extract_dir=os.path.join(folder, 'taxdmp'))
    taxa = sort_taxa_names(folder)
    list_index = get_indeces(folder, taxa)
    nodes_df = get_nodes_file(folder, taxa)
    add_dup_to_taxa(folder, taxa)

