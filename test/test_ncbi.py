import unittest
from unittest.mock import patch
import os 
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import ncbi_tax

class TestGetTaxa(unittest.TestCase): 
    @patch('ncbi_tax.folder', 'test/data')

    def test_load_taxa_file(self): 
        folder = 'test/data'
        taxa = pd.read_csv(os.path.join(folder, 'taxa_names_sorted.tsv'), sep='\t')
        list_index = ncbi_tax.get_indeces(taxa)
        self.assertEqual(27, len(list_index))

        #taxa2, list_index2 = ncbi_tax.get_taxa()
        #self.assertDictEqual(list_index, list_index2)

if __name__=="__main__": 
    unittest.main()