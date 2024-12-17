import pandas as pd
import ncbi_tax
from datetime import datetime
from tqdm import tqdm
import os
import utils

def search_nodes(tax_id:int, nodes_df:pd.DataFrame, mode:str) -> pd.DataFrame: 
    current_tax = tax_id

    lineage = []
    while True: 
        try: 
            node = nodes_df['name_txt'][current_tax]+':'+str(current_tax)
        except KeyError: 
            print('WARNING: Taxon ID '+str(current_tax)+' was not found in the NCBI taxonomy database as stored in the nodes.tsv file.')
            break

        if mode == 'full': 
            lineage.append(node)
        elif mode == 'reduced' and nodes_df['rank'][current_tax] not in ['no rank', 'clade']: 
            lineage.append(node)
        elif mode == 'minimal' and nodes_df.at['rank'][current_tax] in minimal_ranks: 
            lineage.append(node)

        if current_tax == 1: 
            break 
        try: 
            current_tax = nodes_df['parent_tax_id'][current_tax]        
        except KeyError: 
            print('WARNING: Taxon ID '+str(current_tax)+' was not found in the NCBI taxonomy database as stored in the nodes.tsv file.')
            break

    lineage.reverse()

    return lineage

def get_lineage(args): 
    global minimal_ranks

    # Setup
    nodes_df = ncbi_tax.get_nodes()
    output_file = args.prefix+'lineage.tsv'
    minimal_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']

    if args.tax_id: 
        tax_ids = args.tax_id
    elif args.tax_id_file: 
        tax_ids = utils.read_tax_id_file(args.tax_id_file)

    if args.redo or not os.path.exists(output_file): 
        with open(output_file, 'w') as w: 
            w.write('tax_id\tlineage\n')
       
    with open(output_file, 'a') as w: 

        now = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        print(f'\n{now}: Retrieving lineages ({args.lineage})....')

        for tax_id in tqdm(tax_ids, disable=not args.quiet): 
            lineage = search_nodes(tax_id, nodes_df, args.lineage)

            line = str(tax_id)+'\t'
            for lin in lineage: 
                line += lin
                line += ';'

            w.write(line+'\n')     
            if not args.quiet: 
                print(line)

        now = str(datetime.now().strftime("%d/%m/%Y %H:%M:%S"))
        print(f'\n{now}: Results were written into {output_file} file.\n')