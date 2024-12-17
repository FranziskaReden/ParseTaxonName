import ncbi_tax
import get_lineage
import pandas as pd
from fuzzywuzzy import fuzz
import time
from tqdm import tqdm
import multiprocessing
import utils

class Query:
    def __init__(self, name):
        self.original = name
        self.name = utils.tidy_name(name.replace('_', ' '))
        self.red_name = None
        self.min_name = None
        self.tax_id = None
        self.name_txt = None
        self.name_class = None
        self.viral = False
        self.strict_score = 0
        self.relaxed_score = 0
        self.time = None

    def print_info(self): 
        print('original', self.original)
        print('name', self.name)
        print('red_name', self.red_name)
        print('min_name', self.min_name)
        print('tax_id', self.tax_id)
        print('name_txt', self.name_txt)
        print('name_class', self.name_class)
        print('strict_score', self.strict_score)
        print('relaxed_score', self.relaxed_score)

    def update(self, row, score=0):
        self.tax_id = row[0]
        self.name_txt = row[1]
        self.name_class = row[2]
        self.relaxed_score = score
        self.strict_score = fuzz.ratio(self.name, self.name_txt)

    def get_score(self, word:str, candidate:str) -> float: 
        """
        Compares a word with the reference word (candidate) found in the NCBI taxonomy database. 
        Calculates the fuzzy ratio between the word and candidate as well as the original name
        self.name, the reduced named self.red_name as well as the minimal name self.min_name. 
        Returns the highest fuzzy ratio score found.
        """

        candidate = candidate.upper()
        word = word.upper()

        scores = [fuzz.ratio(candidate, self.name.upper()), 0, 0, fuzz.ratio(candidate, word)]
        if self.red_name: 
            scores[1] = fuzz.ratio(candidate, self.red_name.upper())
        if self.min_name: 
           scores[2] = fuzz.ratio(candidate, self.min_name.upper())

        return scores

    def reduce_name(self):
        """
        Tidy up the name (self.name) by omitting some words and replacing others. 
        E.g. 'cf' -> 'cf.', removing all words containing a number, etc. 
        Sets self.red_name (tidied up name) and self.min_name (no numbers) if 
        they are different from self.name.
        """

        # Tidy up name
        original_name = self.name
        original_name = utils.tidy_name2(original_name)

        # Check if sequence is viral
        if any(substring in original_name.upper() for substring in ['VIRAL', 'VIRUS', 'VIRIDAE', 'PHAGE', 'BACTERIOPHAGE']): 
            self.viral = True

        # Return if new name is smaller than 3 (we will not find anything)....
        if len(original_name) < 3:
            return
        
        # Split name
        name_sep = original_name.split(' ')

        # Find words that do not contain useful information        
        name_sep, rm1, rm2 = utils.find_trash_words(name_sep)
        
        # Remove useless words as found in rm1
        for rm in rm1: 
            name_sep.remove(rm)

        new_name = " ".join(name_sep)
        new_name = utils.tidy_name(new_name)

        # Check if the new name is long enough
        if len(new_name)>=3:
            self.red_name = new_name

        else: 
            self.red_name = self.name

        new_name2 = new_name
        for rm in rm2: 
            new_name2 = new_name2.replace(rm, '')
        new_name2 = utils.tidy_name(new_name2)

        # Check if second new name is fferent and long enough
        if len(new_name2) >= 3 and new_name2 != self.red_name: 
            self.min_name = new_name2

class TaxonomySearcher:
    taxa_df = None  # Class-level variable
    list_index = None
    taxa_name_dict = None
    limit = 95

    @classmethod
    def initialize(cls, taxa_df, list_index, taxa_name_dict, limit):
        cls.taxa_df = taxa_df
        cls.list_index = list_index
        cls.taxa_name_dict = taxa_name_dict
        cls.limit = limit

    def __init__(self, name):
        self.name = name

    def get_subset(self, letter):
        if letter in self.list_index:
            limits = self.list_index[letter]
        else:
            limits = self.list_index['_']
        return self.taxa_df.iloc[limits[0]:limits[1]]

    def search_exact(self, query):
        if query.name in self.taxa_name_dict:
            idx = self.taxa_name_dict[query.name]
            query.update(self.taxa_df.iloc[idx].to_numpy())

    def search_approximate(self, query, subset, word):
        matching_indices = subset[subset['name_txt'].str.contains(word, case=False, na=False, regex=False)].index
        best_scores = [0, 0, 0, 0]
        best_candidates = [None, None, None, None]

        for idx in matching_indices:
            candidate = subset.at[idx, 'name_txt']
            if query.viral and not any(v in candidate.upper() for v in ['VIRAL', 'VIRUS', 'VIRIDAE', 'PHAGE', 'BACTERIOPHAGE']):
                continue
            scores = query.get_score(word, candidate)
            for i in range(4):
                if scores[i] > self.limit and scores[i] > best_scores[i]:
                    best_scores[i] = scores[i]
                    best_candidates[i] = idx

        for i in range(4):
            if best_candidates[i] is not None:
                query.update(self.taxa_df.iloc[best_candidates[i]].to_numpy(), best_scores[i])
                break

def start_search(q:Query, searcher:TaxonomySearcher, mode:str): 

    if mode == 'strict': 
        searcher.search_exact(q)
        return
    
    first_letter = q.name[0].upper()
    subset = searcher.get_subset(first_letter)
    
    #relaxed search
    searcher.search_approximate(q, subset, q.name)
    if q.tax_id or mode == 'relaxed': 
        return 

    # lenient search
    # Reduce the name to its most important parts
    q.reduce_name()

    # If we have a reduced name, search with it
    if q.red_name != q.name: 
        if first_letter != q.red_name[0].upper():
            first_letter = q.red_name[0].upper()
            subset = searcher.get_subset(first_letter)
        searcher.search_approximate(q, subset, q.red_name)
    if q.tax_id: 
        return 
    
    # search with the minimal name, if present
    if q.min_name: 
        if first_letter != q.min_name[0].upper():
            first_letter = q.min_name[0].upper()
            subset = searcher.get_subset(first_letter)
        searcher.search_approximate(q, subset, q.min_name)
    if q.tax_id: 
        return 
    
    # Remove words one after another...
    q.min_name = utils.shave_name(q.red_name)
    if q.min_name: 
        if first_letter != q.min_name[0].upper():
            first_letter = q.min_name[0].upper()
            subset = searcher.get_subset(first_letter)

        while q.min_name: 
            searcher.search_approximate(q, subset, q.min_name)
            if q.tax_id: 
                break
            q.min_name = utils.shave_name(q.min_name)

def process_name(args):
    name, mode, searcher = args
    query = Query(name)
    start = time.time()
    start_search(query, searcher, mode)
    query.time = round(time.time() - start, 5)

    if query.tax_id:
        result = (f"{query.original}\t{query.tax_id}\t{query.name_txt}\t{query.name_class}\t"
                  f"{query.strict_score}\t{query.relaxed_score}\t{query.red_name}\t"
                  f"{query.min_name}\t{query.time}")
        return query.tax_id, result
    else:
        return None, query.original

def get_taxids(args):
    taxa_df, list_index = ncbi_tax.get_taxa()
    taxa_name_dict = dict(zip(taxa_df['name_txt'].values, taxa_df.index))

    TaxonomySearcher.initialize(taxa_df, list_index, taxa_name_dict, args.score)
    searcher = TaxonomySearcher('ncbi')

    output_files = args.prefix + "tax_ids.tsv", args.prefix + 'tax_ids_failed.txt'
    processed_names, failed_names = utils.load_checkpoint(output_files, args.redo)

    # Names to process
    if args.taxon_name:
        names = args.taxon_name
    elif args.ali_file:
        print(f"Reading in alignment file {args.ali_file}...")
        names = utils.read_ali_file(args.ali_file)
    elif args.name_file:
        print(f"Reading in name file {args.name_file}...")
        names = utils.read_name_file(args.name_file)

    unique_names = set(names)
    print(f"Loaded {len(names)} names. Of those {len(unique_names)} are unique.")

    # Exclude already processed names
    if len(processed_names) > 0: 
        names_to_process = [name for name in unique_names if name not in processed_names and name not in failed_names]
    else: 
        names_to_process = unique_names

    if len(names_to_process) == 0: 
        print(f'0 new names to process were found. Matched and failed names can be found in files \
              {output_files[0]} and {output_files[1]} repectivly. \
              \nUse the --redo flag should you wish to rerun the analysis, which will overwrite the \
              results file.')
        return 

    results, failed, tax_ids = [], [], []
    for name in tqdm(names_to_process, disable=not args.quiet):
        result = process_name((name, 'strict', searcher))
        if result[0] is None:
            failed.append(result[1])
        else:
            results.append(result[1])
            tax_ids.append(int(result[0]))

    # Final checkpoint save
    utils.write_checkpoint(output_files, results, failed, len(results), mode = False, quiet=args.quiet)
    searcher.taxa_name_dict.clear()

    # Approximate or lenient search
    processed_count = 0
    if args.mode != 'strict' and len(failed) > 0:

        failed2 = []

        # Prepare multiprocessing
        if args.cores in ['AUTO', 'auto']: 
            num_processes = multiprocessing.cpu_count()
        else: 
            num_processes = int(args.cores)

        pool = multiprocessing.Pool(num_processes)

        # Prepare arguments for parallel processing
        args_list = [(name, args.mode, searcher) for name in failed]
        if not args.quiet:
            print("name\ttax_id\tname_txt\tname_class\tstrict_score\t"
            "relaxed_score\treduced_name\tmin_name\ttime(s)")
        with pool:
            for result in tqdm(pool.imap(process_name, args_list), total=len(failed), disable=not args.quiet):
                processed_count += 1

                if result[0] is None:
                    failed2.append(result[1])  # Collect failed2 names
                else:
                    results.append(result[1])
                    tax_ids.append(int(result[0]))

                    if not args.quiet: 
                        print(result[1])

                # Periodically save checkpoint
                if processed_count % 500 == 0:
                    utils.write_checkpoint(output_files, results, failed2, processed_count, mode=True, quiet=args.quiet)

        # Final checkpoint save
        utils.write_checkpoint(output_files, results, failed2, processed_count, mode=True)

    print(f"Matched names written to {output_files[0]}. Failed names to {output_files[1]}.")
    if tax_ids:
        args.tax_id = tax_ids
        get_lineage.get_lineage(args)