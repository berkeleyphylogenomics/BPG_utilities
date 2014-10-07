'''
    Okay, so none of us know django all that great, and we don't really have test cases.
    Our schema breaks the db creation jazz, so using the django stuff is less trivial
    than I have time to figure out right now.
'''

import unittest
from annotation import TreeAnnotation, GOAnnotation, consensus, summarize

class TestTreeAnnotation(unittest.TestCase):
    def test_echo(self):
        a = TreeAnnotation([1,2,3])
        self.assertEqual(a.annotation, consensus([a]).annotation)

    def test_shorten(self):
        a = TreeAnnotation([1,2,3])
        b = TreeAnnotation([1,4,3])
        self.assertEqual([1], consensus([a, b]).annotation)

    def test_summarize(self):
        a = TreeAnnotation([1,2,3])
        b = TreeAnnotation([1,2,3])
        c = TreeAnnotation([0,2,3])
        self.assertEqual(2, len(summarize([a, b, c])))

    def test_ec(self):
        a = TreeAnnotation([
            { 'name' : 'Hydrolases.', 'value': '3' },
            { 'name' : 'Acting on ester bonds.', 'value': '1' },
            { 'name' : "Endoribonucleases producing 5'-phosphomonoesters.", 'value': '26' },
            { 'name' : 'Ribonuclease E.', 'value': '12' },
            ])
        b = TreeAnnotation([
            { 'name' : 'Hydrolases.', 'value': '3' },
            { 'name' : 'Acting on ester bonds.', 'value': '1' },
            { 'name' : "Endoribonucleases producing 5'-phosphomonoesters.", 'value': '26' },
            { 'name' : 'Ribonuclease (poly-(U)-specific).', 'value': '9' },
            ])
        self.assertEqual(consensus([a, b]).annotation,
            [{ 'name' : 'Hydrolases.', 'value': '3' },
            { 'name' : 'Acting on ester bonds.', 'value': '1' },
            { 'name' : "Endoribonucleases producing 5'-phosphomonoesters.", 'value': '26' }])

        
class TestGOAnnotation(unittest.TestCase):
    def test_accum(self):
        a = GOAnnotation({
            'accession': 'GO:0019143',
            'evidence_priority': 16
        })
        b = GOAnnotation({
            'accession': 'GO:0008690',
            'evidence_priority': 5
        })
        self.assertEqual([ a.annotation['evidence_priority'] for a in summarize([a, b]) ], [5, 16])

    def test_merge(self):
        a = GOAnnotation({
            'accession': 'GO:0008690',
            'evidence_priority': 16
        })
        b = GOAnnotation({
            'accession': 'GO:0008690',
            'evidence_priority': 5
        })
        self.assertEqual([ a.annotation['evidence_priority'] for a in summarize([a, b]) ], [5])

if __name__ == '__main__':
    unittest.main()
