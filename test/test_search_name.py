import unittest
from unittest.mock import patch
import os 
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import search_name as sn
import ncbi_tax, utils

class TestGetTaxa(unittest.TestCase): 
    @patch('ncbi_tax.folder', f'{os.path.join(os.environ.get("HOME"), ".ncbi_tax")}')

    def test_query_class(self): 
        
        name = '345_uncultured_eukaryote_SLV_3GJ1_11_KT072099'
        q = sn.Query(name)

        self.assertEqual(q.original, name)
        self.assertEqual(q.name, name.replace('_', ' '))

        q.reduce_name()
        self.assertEqual(q.red_name, 'uncultured_eukaryote_SLV_3GJ1_11_KT072099'.replace('_', ' '))
        self.assertEqual(q.min_name, 'uncultured_eukaryote_SLV'.replace('_', ' '))

        scores = q.get_score('uncultured eukaryote', 'Uncultured Eukaryote')
        self.assertEqual([62, 66, 91, 100], scores)

    def test_initialization(self): 

        taxa_df, list_index = ncbi_tax.get_taxa()
        taxa_name_dict = dict(zip(taxa_df['name_txt'].values, taxa_df.index))
        score = 95
        mode = 'lenient'

        sn.TaxonomySearcher.initialize(taxa_df, list_index, taxa_name_dict, score)
        searcher = sn.TaxonomySearcher('ncbi')

        self.assertEqual(score, searcher.limit)
        self.assertDictEqual(list_index, searcher.list_index)

        # Example for exact search
        name = 'Homo_sapiens '
        tax_id = 9606 

        q = sn.Query(name)
        self.assertEqual(name, q.original)
        self.assertEqual(name.replace('_', ' ').strip(), q.name)

        # Test exact search function
        searcher.search_exact(q)
        self.assertEqual(tax_id, q.tax_id)

        # Test exact search in embedded in start_search function
        q = sn.Query(name)
        sn.start_search(q, searcher, 'strict')
        self.assertEqual(tax_id, q.tax_id)

        # Test approximate search
        name = 'Mus musculu'

        q = sn.Query('345_uncultured_eukaryote_SLV_3GJ1_11_KT072099')
        sn.start_search(q, searcher, mode)

        subset = searcher.get_subset(q.name[0].upper())
        searcher.search_approximate(q, subset, 'uncultured eukaryote')
        

if __name__=="__main__": 
    unittest.main()
