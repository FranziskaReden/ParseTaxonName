import unittest
from unittest.mock import patch
import os 
import sys
import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
import ncbi_tax, get_lineage

class TestFerLineage(unittest.TestCase): 

    def test_search_nodes(self):

        folder = f'{os.path.join(os.environ.get("HOME"), ".ncbi_tax")}'

        nodes_df = ncbi_tax.get_nodes(folder)

        reduced_ranks = ['superkingdom', 'genus', 'species', 'order', 'family',
       'subspecies', 'subfamily', 'strain', 'serogroup', 'tribe',
       'phylum', 'class', 'species group', 'forma', 'subphylum',
       'suborder', 'subclass', 'varietas', 'kingdom', 'forma specialis',
       'isolate', 'superfamily', 'infraorder', 'infraclass', 'superorder',
       'subgenus', 'superclass', 'parvorder', 'serotype',
       'species subgroup', 'subcohort', 'cohort', 'subtribe', 'section',
       'series', 'subkingdom', 'superphylum', 'subsection']
        minimal_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']

        tax_id = 9606

        mode = 'full'
        lineage = get_lineage.search_nodes(tax_id, nodes_df, [reduced_ranks, minimal_ranks], mode)
        line = ''
        for lin in lineage: 
            line += lin
            line += ';'
        self.assertEqual('root:1;cellular organisms:131567;Eukaryota:2759;Opisthokonta:33154;Metazoa:33208;Eumetazoa:6072;Bilateria:33213;Deuterostomia:33511;Chordata:7711;Craniata:89593;Vertebrata:7742;Gnathostomata:7776;Teleostomi:117570;Euteleostomi:117571;Sarcopterygii:8287;Dipnotetrapodomorpha:1338369;Tetrapoda:32523;Amniota:32524;Mammalia:40674;Theria:32525;Eutheria:9347;Boreoeutheria:1437010;Euarchontoglires:314146;Primates:9443;Haplorrhini:376913;Simiiformes:314293;Catarrhini:9526;Hominoidea:314295;Hominidae:9604;Homininae:207598;Homo:9605;Homo sapiens:9606;', line)

        mode='reduced'
        lineage = get_lineage.search_nodes(tax_id, nodes_df, [reduced_ranks, minimal_ranks], mode)
        line = ''
        for lin in lineage: 
            line += lin
            line += ';'
        self.assertEqual('root:1;Eukaryota:2759;Metazoa:33208;Chordata:7711;Craniata:89593;Sarcopterygii:8287;Mammalia:40674;Euarchontoglires:314146;Primates:9443;Haplorrhini:376913;Simiiformes:314293;Catarrhini:9526;Hominoidea:314295;Hominidae:9604;Homininae:207598;Homo:9605;Homo sapiens:9606;', line)

        mode='minimal'
        lineage = get_lineage.search_nodes(tax_id, nodes_df, [reduced_ranks, minimal_ranks], mode)
        line = ''
        for lin in lineage: 
            line += lin
            line += ';'
        self.assertEqual('root:1;Eukaryota:2759;Metazoa:33208;Chordata:7711;Mammalia:40674;Primates:9443;Hominidae:9604;Homo:9605;Homo sapiens:9606;', line)
        
