import unittest
import os 
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import utils

class TestUtils(unittest.TestCase): 

    def test_read_line(self):
        self.assertEqual(utils.read_line('144551	|	 Krishnamurthy 11-00121	|	 Krishnamurthy 11-00121 <holotype>	|	type material	|'), \
                         ['144551', 'Krishnamurthy 11-00121', ' Krishnamurthy 11-00121 <holotype>', 'type material']) 
        self.assertEqual(utils.read_line('144509	|	Fusarium sp. BBA 65925	|		|	scientific name	|'), \
                         ['144509', 'Fusarium sp. BBA 65925', '', 'scientific name']) 
        
    def test_tidy_name(self): 
        self.assertEqual(utils.tidy_name(' Homo     sapiens '), 'Homo sapiens')

    def test_tidy_name2(self): 
        self.assertEqual(utils.tidy_name2('Homo_sapiens'), 'Homo sapiens')
        self.assertEqual(utils.tidy_name2('[Homo]_sapiens'), ' Homo  sapiens')
        self.assertEqual(utils.tidy_name2('Homo.sapiens.'), 'Homo sapiens ')
        self.assertEqual(utils.tidy_name2('Homo-sapiens'), 'Homo sapiens')
        self.assertEqual(utils.tidy_name2('Virus_ex_Homo_sapiens'), 'Virus')

    def test_has_number(self): 
        self.assertEqual(utils.has_number('CMW10125'), True)

    def test_find_trash_words(self): 
        self.assertEqual(utils.find_trash_words(['Mcclungia', 'sp', 'GENOME', '135efg','cymo', 'ME11_159']), \
                         (['Mcclungia', 'sp.', 'GENOME', '135efg', 'cymo', 'ME11_159'],['GENOME'],['135efg', 'ME11_159']))
        
    def test_read_name_file(self): 
        self.assertEqual(utils.read_name_file('test/data/names_list.txt'), (['Homo sapiens', 'Mus musculus']))
    
    def test_checkpoint(self): 

        file_path = ['test/data/test_checkpoint_tax_ids.txt', 'test/data/test_checkpoint_tax_ids_failed.txt']        
        results = ['Homo sapiens\t9606\tHomo sapiens\tscientific name\t100\t0\tNone\tNone\t0.00016', \
        'Drosophila sp\t7242\tDrosophila sp.\tincludes\t96\t96\tNone\tNone\t0.03541', \
            'Mus muskulus\t10090\tMus musculus\tscientific name\t92\t92\tMus muskulus\tMus\t0.09495']
        failed = ['Homo sapiens', 'Drosophila sp', 'Mus muskulus']

        self.assertEqual(utils.load_checkpoint(file_path, True), (set(), set()))
        utils.write_checkpoint(file_path, results.copy(), failed.copy(), 3, mode = False)
        self.assertEqual(utils.load_checkpoint(file_path, False), ({'Mus muskulus', 'Homo sapiens', 'Drosophila sp'}, set()))
        utils.load_checkpoint(file_path, True)
        utils.write_checkpoint(file_path, results.copy(), failed.copy(), 3, mode = True)
        self.assertEqual(utils.load_checkpoint(file_path, False), ({'Mus muskulus', 'Homo sapiens', 'Drosophila sp'}, set(failed)))

    def test_read_taxid_file(self): 
        self.assertEqual(utils.read_tax_id_file('test/data/tax_id_list.txt'), [9606, 10090])

    def test_shave_name(self): 
        word = 'mu uncultured eukaryote SLV 3GJ1 11 KT072099'
        word_less = 'mu uncultured eukaryote SLV 3GJ1 11'
        self.assertEqual(utils.shave_name(word), word_less)

        word = 'mu uncultured eukaryote SLV 3GJ1 11'
        word_less = 'mu uncultured eukaryote SLV 3GJ1'
        self.assertEqual(utils.shave_name(word), word_less)

        word = 'mu uncultured eukaryote SLV 3GJ1'
        word_less = 'mu uncultured eukaryote SLV'
        self.assertEqual(utils.shave_name(word), word_less)

        word = 'mu uncultured eukaryote SLV'
        word_less = 'mu uncultured eukaryote'
        self.assertEqual(utils.shave_name(word), word_less)
        
        word = 'mu uncultured'
        word_less = None
        self.assertEqual(utils.shave_name(word), word_less)
        

if __name__=="__main__": 
    unittest.main()

