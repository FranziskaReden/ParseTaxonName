import os

def shave_name(word:str) -> str | None:
    '''removes last word from string. Returns reduced 
    name or None, if new string too small (<=3) or 
    if it was the last word in string.'''

    new_name = None
    split_word = word.split(' ')

    if word and len(split_word) > 1:
        split_word.remove(split_word[-1])
        new_name = rejoin_name(split_word)
        if len(new_name) >= 3:
            return new_name

    return None

def read_line(line: str) -> list:
    '''Reads in line seperated with "|". Returns a list of 
    strings containing the cleaned up entries of the "|" 
    seperated file'''

    line = line.strip('\n')
    line = line.replace('\t', '')
    line = line.split('|')
    line[1] = line[1].replace('\"', '')
    line[1] = line[1].replace('\'', '')
    if line[-1] == '':
        line = line[:-1]
    if line[1][0] == ' ':
        line[1] = line[1][1:]

    return line

def tidy_name(word:str) -> str:
    '''Removes whitespaces from input string word. 
    Returns tidied up word.'''
    while '  ' in word:
        word = word.replace("  ", " ")
    return word.strip()

def tidy_name2(original_name:str) -> str:
    '''Removes whitespaces from input string word. 
    Returns tidied up word.'''

    if '_ex_' in original_name:
        original_name = original_name.split('_ex_')[0]

    original_name = original_name.replace('_', ' ')
    if original_name[0] != '[':
        original_name = original_name.replace('[', '(')
        original_name = original_name.replace(']', ')')

    if original_name.count('.') > 1:
        original_name = original_name.replace('.', ' ')
    if original_name.count('-') > 1:
        original_name = original_name.replace('-', ' ')

    return original_name

def has_number(word:str) -> bool:
    '''Function to check if any character in a given string (word) 
    is numeric. Retruns True or False accordingly'''
    return any(str.isnumeric(c) for c in word)

def rejoin_name(name_sep:list) -> str:

    new_name = " ".join(name_sep)
    return tidy_name(new_name)

def find_trash_words(name_sep:list) -> list:
    '''Check each word in the input list of strings name_sep. 
    Find words to remove from string that are not needed.
    Return the name, and two sets of list containing the words 
    to be removed.'''

    rm1 = []
    rm2 = []
    tmp = 0

    for i, value in enumerate(name_sep):
    # Remove all words containing numbers until encountering first word that does not.
        if tmp == 0 and has_number(name_sep[i]) is True:
            rm1.append(name_sep[i])
        else:
            tmp = 1

        if name_sep[i] in ['sp', 'cf', 'pv', 'aff', 'var']:
            name_sep[i] = name_sep[i] + '.'

        elif name_sep[i].upper() in ['GENOME', 'REVERSED']:
            rm1.append(name_sep[i])

        if has_number(name_sep[i]) is True or len(name_sep[i])<3:
            if name_sep[i] not in rm1:
                rm2.append(name_sep[i])

    return name_sep, rm1, rm2

def read_ali_file(file_name:str) -> list:
    '''To be implemented....'''

    names = []

    return names

def read_name_file(file_name:str) -> list:
    '''Reads in file file_name. Returns list containing the lines of the file.'''

    with open(file_name, encoding='utf-8') as t:
        names = t.readlines()

    for i, value in enumerate(names):
        names[i] = names[i].strip()

    return names

def write_checkpoint(file_path: str, results: list, failed: list, processed_count: int,
                     mode:bool = False, quiet:bool=False) -> None:
    """Write progress to a checkpoint file.
    
    Parameters
    ----------
    file_path : str
        Path to the file to write the ceckpoint to
    results : list
        list of sucessful queries
    failed : list
        list of failed queries
    processed_count : int 
        Number of queries processed
    mode : bool
        Mode declares whether to clear both results and 
        failed lists after writing out the resuluts or 
        just results. default=False
    quiet : bool
        States whether to print process to the screen 
        Default=False
    """

    checkpoint_path = file_path[0]
    failed_path = file_path[1]

    with open(checkpoint_path, "a", encoding='utf-8') as w:
        for result in results:
            if result:
                w.write(result + "\n")

    if mode is True:
        with open(failed_path, "a", encoding='utf-8') as f:
            for name in failed:
                f.write(name + "\n")

    if quiet is False:
        print(f"Checkpoint saved: {processed_count} names processed.")

    # Clear the lists to free memory
    results.clear()
    if mode:
        failed.clear()

def load_checkpoint(file_path, redo):
    """Load processed names from checkpoint files."""

    processed_names = set()
    failed_names = set()

    checkpoint_path = file_path[0]
    failed_path = file_path[1]

    if os.path.exists(checkpoint_path) and not redo:
        with open(checkpoint_path, "r", encoding='utf-8') as r:
            next(r)
            for line in r:
                processed_names.add(line.split("\t")[0])
    else:
        with open(checkpoint_path, "w", encoding='utf-8') as w:
            w.write("name\ttax_id\tname_txt\tname_class\tstrict_score\t"
            "relaxed_score\treduced_name\ttmp_name\tmin_name\ttime(s)\n")

    if os.path.exists(failed_path) and not redo:
        with open(failed_path, "r", encoding='utf-8') as f:
            for line in f:
                failed_names.add(line.strip())

    elif os.path.exists(failed_path):
        os.remove(failed_path)

    if len(processed_names) > 0 or len(failed_names) > 0:
        print(f"\nLoaded {len(processed_names)} processed \
              names and {len(failed_names)} failed names from checkpoint.")

    return processed_names, failed_names

def read_tax_id_file(file_name:str) -> list:
    ''''Read in file containing tax_ids.
    Returns list of tax_ids (as integers).'''

    with open(file_name, encoding='utf-8') as t:
        tax_ids = t.readlines()

    for i, value in enumerate(tax_ids):
        tax_ids[i] = int(tax_ids[i].strip())

    return tax_ids
