import argparse
import pathlib
import os
import sys

import search_name
import get_lineage
import ncbi_tax

def main():
    '''Script to retrieve taxon ID according \
    to te NCBI taxonomy given a taxon name and/or to \
    retrieve the lineage of said/a taxon ID.
    '''

    parser = argparse.ArgumentParser(description='**Script to retrieve taxon ID according \
                            to te NCBI taxonomy given a taxon name and/or to \
                            retrieve the lineage of said/a taxon ID.**')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-n', '--taxon_name', type=str, action='store', nargs='+',
                    help='Taxon name for which to find the taxon ID')
    group.add_argument('-id', '--tax_id', type=int, action='store', nargs='+',
                        help='Taxon ID for which to return the lineage')
    group.add_argument('-a', '--ali_file', type=pathlib.Path, action='store',
                        help='Alignment file. Will find taxon IDs of the sequence \
                            names given in the file.')
    group.add_argument('-nf', '--name_file', type=pathlib.Path, action='store',
                        help='File containing a list of taxon names for which to find \
                            the NCBI taxon ID. Each name has to be stated in a new line.')
    group.add_argument('-idf', '--tax_id_file', type=pathlib.Path, action='store',
                        help='File containing a list of taxon IDs for which to extract the \
                            lineage. Each taxon ID has to be stated in a new line.')
    parser.add_argument('-l', '--lineage', type=str, choices=['full','reduced','minimal'],
                        action='store', default='full',
                        help='States whether to return the full, reduced (only unique ranks), \
                            or minimal lineage. Default is full.')
    parser.add_argument('--mode', type=str, choices=['strict','relaxed','lenient'],
                        default='strict', help='States how strict the search for taxon ID \
                            given a taxon name should be. Strict mode will only return \
                            perfect matches. Relaxed mode will allow for mismatches. \
                            Lenient mode will potentially modify the given taxon name \
                            (e.g., remove accession number from the name) to only retain the \
                            informative parts. Default is strict mode.')
    parser.add_argument('--score', type=float, default=95,
                        help='If mode is relaxed or lenient, score states the minimal matching \
                            score between the given name and a matching name found in the NCBI \
                            taxonomy. Can range from 60-100. Default is 95.')
    parser.add_argument('--prefix', type=str, default='', help='Prefix for the output files.')
    parser.add_argument('-r', '--redo', default=False, action = 'store_true',
                        help='Will redo analysis. Be aware that output files will be \
                            overwritten.')
    parser.add_argument('-q', '--quiet', default=False, action = 'store_true',
                        help='Quiet mode, will output minimal information on your screen.')
    parser.add_argument('--cores', default=1, action='store', type=int,
                        help='For parallellized searching. Declares the number of cores to us. \
                            Default is 1.')
    parser.add_argument('-db', default=str(os.path.join(os.environ.get("HOME"), ".ncbi_tax")),
                        action='store',
                        help='Path to and name of the folder in which to write/find the \
                            NCBI taxonomy database. Default is the older ".ncbi_tax" in \
                            the home directory of the current user.')
    parser.add_argument('--update', default=False, action='store_true',
                        help='Will update NCBI taxonomy database (will be downloaded).')
    args = parser.parse_args()

    if args.score > 100 or args.score < 60:
        print('Please choose a score value between 60 and 100! \
              I do not recommend going below 90%.\n \
              Please, also note that a value of 100 in relaxed \
              mode equals strict mode (which is more efficient).\n')
        sys.exit(2)

    print(parser.description, '\n')

    if args.prefix != '':
        args.prefix = args.prefix+'_'
    elif args.ali_file:
        args.prefix = str(args.ali_file) + '_'
    elif args.name_file:
        args.prefix = str(args.name_file) + '_'
    elif args.tax_id_file:
        args.prefix = str(args.tax_id_file) + '_'

    if args.update is True:
        ncbi_tax.update_db(args.db)

    if args.taxon_name or args.ali_file or args.name_file:
        search_name.get_taxids(args)

    elif os.path.isfile(args.prefix+'lineage.tsv') and args.redo is False:
        print('Output file '+args.prefix+'_lineage.tsv detetced. Nothing to do... \
              \nUse the --redo flag should you wish to rerun the analysis, \
              which will overwrite the results file.\n')

    else:
        get_lineage.get_lineage(args)

if __name__ == "__main__":
    main()
