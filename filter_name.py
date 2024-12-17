import pandas as pd
import sys
from difflib import SequenceMatcher
import copy
import os.path
import shutil
import string
import urllib
import urllib.request
from socket import timeout
from Bio import SeqIO

def reduce_name(name:str) -> list: 
    '''
    Function to 'tidy up' the sequence label as found in an alignments file, which is given as input. 

    The name is seperated into its words. All words containing numbers (possible accession numbers) as 
    found at the beginning of the name are removed from the name. Abbrevations such as 'sp.' that are missing 
    a point at the end, are replaced (e.g., sp -> sp.). 

    INPUT
    --------
    name:str  
        The sequence name to be 'tidied up'. 

    RETURNS
    --------
    Output is a list:
        - At position 0 of the list the tidied up name can be found. The new name can also be None if the name is too short 
        (less than 3 characters as it contains too little information) or if the name does not contain any 
        useful information (e.g., only numbers). 
        - At position 1, an alternative 'tidied up' name can be found. It alternatively 
        removes all words containing digits as well as words of length<=2. 
        - Position 2 of the returned list contains a boolean variable indicating wether there are hints towards the sequence being viral (if yes -> True). 
    '''

    viral = False

    # If the string '_ex_' (standing for extracted from) can be found in the name, only consider the name after '_ex_'
    if '_ex_' in name: 
        name = name.split('_ex_')[1]
    name = name.replace('_', ' ')

    # Check if the sequence name contains any indication for the sequence to be viral
    if any(substring in name.upper() for substring in ['VIRAL', 'VIRUS', 'VIRIDAE']): 
        viral = True

    # If the name is shorter than 3 characters, return None (too little information to search)
    if len(name)<3:
        return None, None, viral
    # Check if points are the delimiter (instead of underscore)
    if name.count('.') > 1: 
        name = name.replace('.', ' ')
    if '-' in name: 
        name = name.replace('-', ' ')
    # Split name into words
    name_sep = name.split(' ')

    # Declare a list that will contain all words that will be removed from the name
    to_remove_1 = []
    to_remove_2 = []
    tmp = 0
    for i in range (len(name_sep)): 
        # Remove all words containing numbers until encountering first word that does not.
        if tmp == 0 and has_numbers(name_sep[i]) is True:
            to_remove_1.append(name_sep[i]) 
        else: 
            tmp = 1  
        if name_sep[i] in ['sp', 'cf', 'pv', 'aff', 'var']: 
            name_sep[i] == name_sep[i]+'.'
        elif name_sep[i].upper() in ['GENOME', 'REVERSED']: 
            to_remove_1.append(name_sep[i])      
        if has_numbers(name_sep[i]) is True or len(name_sep[i])<3:
            if name_sep[i] not in to_remove_1: 
                to_remove_2.append(name_sep[i]) 

    for rm in to_remove_1: 
        name_sep.remove(rm)

    # Put all words into a single string again    
    new_str = ''
    for i in range (len(name_sep)): 
        if i < (len(name_sep)-1): 
            new_str+=name_sep[i]+' '
        else: 
            new_str+=name_sep[i] 

    new_str_2 = copy.copy(new_str)
    # Create second name that additionally contains no words with numbers.
    for rm in to_remove_2: 
        new_str_2 = new_str_2.replace(rm, '')
    
    while '  ' in new_str: 
        new_str = new_str.replace('  ', ' ')
    while len(new_str)>0 and new_str[0]==' ': 
        new_str = new_str[1:]
    while len(new_str)>0 and new_str[-1]==' ': 
        new_str = new_str[:-1]

    while '  ' in new_str_2: 
        new_str_2 = new_str_2.replace('  ', ' ')
    while len(new_str_2)>0 and new_str_2[0]==' ': 
        new_str_2 = new_str_2[1:]
    while len(new_str_2)>0 and new_str_2[-1]==' ': 
        new_str_2 = new_str_2[:-1]

    # If the new name is too short or if the original name does not hold any relevent information, return None
    if len(new_str)<3 or new_str == '': 
        return None, None, viral
    
    # Else, return the new name(s).
    if len(new_str_2)<3 or new_str_2 == '': 
        return new_str, None, viral
    else: 
        return new_str, new_str_2, viral