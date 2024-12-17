import os

def shave_name(word:str): 

    new_name = None
    split_word = word.split(' ')
        
    if word and len(split_word) > 1: 
        new_name = word.replace(' '+split_word[-1], '')
        if len(new_name) >= 3: 
            return new_name

    return None

def read_line(line: str): 

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

def tidy_name(word:str): 
    while '  ' in word: 
        word = word.replace("  ", " ")
    return word.strip()
    
def tidy_name2(original_name:str): 

    if '_ex_' in original_name:
        original_name = original_name.split('_ex_')[0]

    original_name = original_name.replace('_', ' ')
    original_name = original_name.replace('[', ' ')
    original_name = original_name.replace(']', ' ')

    if original_name.count('.') > 1: 
        original_name = original_name.replace('.', ' ')
    if '-' in original_name: 
        original_name = original_name.replace('-', ' ')

    return original_name

def has_number(word:str): 
    return any(str.isnumeric(c) for c in word)

def find_trash_words(name_sep:str):

    rm1 = []
    rm2 = []
    tmp = 0

    for i in range (len(name_sep)): 
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

    names = []

    return names

def read_name_file(file_name:str) -> list:

    with open(file_name) as t: 
        names = t.readlines()

    for i in range (len(names)): 
        names[i] = names[i].strip()

    return names

def write_checkpoint(file_path, results, failed, processed_count, mode = False, quiet=False):
    """Write progress to a checkpoint file."""
    checkpoint_path = file_path[0]
    failed_path = file_path[1]

    with open(checkpoint_path, "a") as w:
        for result in results:
            if result:
                w.write(result + "\n")

    if mode is True: 
        with open(failed_path, "a") as f:
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
        with open(checkpoint_path, "r") as r:
            next(r)
            for line in r:
                processed_names.add(line.split("\t")[0])  # Assuming first column is the name
    else: 
        with open(checkpoint_path, "w") as w:
            w.write("name\ttax_id\tname_txt\tname_class\tstrict_score\t"
            "relaxed_score\treduced_name\tmin_name\ttime(s)\n")

    if os.path.exists(failed_path) and not redo:
        with open(failed_path, "r") as f:
            for line in f:
                failed_names.add(line.strip())

    elif os.path.exists(failed_path): 
        os.remove(failed_path)

    if len(processed_names) > 0 or len(failed_names) > 0: 
        print(f"\nLoaded {len(processed_names)} processed names and {len(failed_names)} failed names from checkpoint.")

    return processed_names, failed_names  

def read_tax_id_file(file_name:str) -> list:

    with open(file_name) as t: 
        tax_ids = t.readlines()
    
    for i in range (len(tax_ids)): 
        tax_ids[i] = int(tax_ids[i].strip())

    return tax_ids