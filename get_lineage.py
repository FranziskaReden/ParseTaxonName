import os
from datetime import datetime
from tqdm import tqdm
import pandas as pd
from collections import Counter

import utils
import ncbi_tax

def search_nodes(tax_id:int, nodes_df:pd.DataFrame, ranks:list, mode:str) -> list:
    '''Function to retrieve the lineage of a given taxon ID. 
    Returns lineage.
    
    Parameters
    ----------
    tax_id : int
        Taxon ID for which to find the lineage.
    nodes_df : pd.DataFrame 
        DataFrame holding the nodes of the NCBI taxonomy database
    ranks : list 
        List of lists: list of reduced and minimal ranks respectively. 
    mode : str 
        States whether to return full, reduced or minimal lineage.
    
    Returns:
    ----------
    lineage : list 
        List of strigs containing the lineage of the given tax_id.
    '''

    reduced_ranks, minimal_ranks = ranks
    current_tax = tax_id

    lineage = []
    while True:
        try:
            node = nodes_df.at[current_tax, 'name_txt']+':'+str(current_tax)
        except KeyError:
            print('WARNING: Taxon ID '+str(current_tax)+' was not found in the NCBI \
                  taxonomy database as stored in the nodes.tsv file.')
            break

        if mode == 'full':
            lineage.append(node)
        elif mode == 'reduced' and nodes_df.at[current_tax, 'rank'] \
            in reduced_ranks or current_tax == 1:
            lineage.append(node)
        elif mode == 'minimal' and nodes_df.at[current_tax, 'rank'] \
            in minimal_ranks or current_tax == 1:
            lineage.append(node)

        if current_tax == 1:
            break
        try:
            current_tax = nodes_df['parent_tax_id'][current_tax]
        except KeyError:
            print('WARNING: Taxon ID '+str(current_tax)+' was not found in the NCBI \
                  taxonomy database as stored in the nodes.tsv file.')
            break

    lineage.reverse()

    return lineage

def get_lineage(args):
    '''Function to retrieve the lineage given the arguments parsed from the 
    command line in the main function. Results are writen into an output file. 
    Returns None.
    
    Parameters
    ----------
    args: argparse.Namespace
        arguments parsed from command line in main function. 
    '''

    # Setup
    nodes_df = ncbi_tax.get_nodes(args.db)
    output_file = args.prefix+'lineage.tsv'
    reduced_ranks = ['domain', 'kingdom', 'subkingdom', 'superphylum', 
                     'subphylum', 'phylum', 'superclass', 'class', 'subclass', 
                     'infraclass', 'cohort', 'subcohort', 'superorder', 'order', 
                     'suborder', 'infraorder', 'parvorder', 'superfamily', 'family', 
                     'subfamily', 'genus', 'subgenus', 'species group', 'species subgroup', 
                     'species', 'subspecies', 'tribe', 'subtribe', 'forma', 'varietas', 
                     'strain', 'section', 'subsection', 'pathogroup', 'subvariety', 
                     'genotype', 'serotype', 'isolate', 'morph', 'series', 
                     'forma specialis', 'serogroup', 'biotype', 'acellular_root', ]
    minimal_ranks = ['species', 'genus', 'family', 'order', 'class',
        'phylum', 'kingdom', 'domain']

    tax_ids = []

    if args.tax_id:
        tax_ids = args.tax_id
    elif args.tax_id_file:
        tax_ids = utils.read_tax_id_file(args.tax_id_file)

    if args.redo or not os.path.exists(output_file):
        with open(output_file, 'w', encoding='utf-8') as w:
            w.write('tax_id\tquantity\tlineage\n')

    with open(output_file, 'a', encoding='utf-8') as w:

        now = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        print(f'\n{now}: Retrieving lineages ({args.lineage})....')
        
        unique_tax_ids = Counter(tax_ids)

        for tax_id, count in tqdm(unique_tax_ids.items(), disable=not args.quiet):
            lineage = search_nodes(tax_id, nodes_df, [reduced_ranks, minimal_ranks], args.lineage)
            line = f'{tax_id}\t{count}\t'
            for lin in lineage:
                line += f'{lin};'
            w.write(line+'\n')

            if not args.quiet:
                print(line)

        now = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        print(f'\n{now}: Results were written into {output_file} file.\n')
