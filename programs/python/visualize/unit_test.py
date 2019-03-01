import unittest
import visualize_nucleotide_counts as vnc


class TestFileReaders(unittest.TestCase):

    def test_BuildDomTblList(self):
        testList = vnc.buildDomTblList("test/output/domtbl/Pf_prophage_LESB58.domtbl")
        self.assertEqual(testList[0], ['TIGR02224', '9879', '1.2e-23', 'recomb_XerC: tyrosine recombinase XerC'])

    def test_BuildHmmStatDict(self):
        testDict = vnc.buildHmmStatDict("")


if __name__ == '__main__':
    unittest.main()
