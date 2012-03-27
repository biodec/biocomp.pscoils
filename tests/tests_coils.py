#import sys
#sys.path.insert('src')
from unittest import TestCase
import StringIO
from biocomp.pscoils import pred_coil
from biocomp.pscoils.utils import readFasta, printCoil
from biocomp.pscoils.Params import Params
from Bio.SeqRecord import SeqRecord

class test_suite(TestCase):
    def test_ferritin(self):
        # ../psCoils.py -f ./ferritin.fasta -l F > ./ferritin.Coils
        weight='un'
        wlambda=0.5
        win=21
        modeprofile=False
        seqs=readFasta("tests/ferritin.fasta")
        modefasta=True
        labels='F'
        params=Params(weight,win)
        output = StringIO.StringIO()
        for name in seqs.keys():
            seqr = SeqRecord(seqs[name], id=name)
            seqr = pred_coil(seqr, params)
            printCoil(seqr,labels,io=output)
        bench_output = open("tests/ferritin.Coils").read()
        self.assertEqual(output.getvalue(),bench_output)        
        
    def test_ceba_bovin(self):
        # ../psCoils.py -f ./ferritin.fasta -l F > ./ferritin.Coils
        weight='un'
        wlambda=0.5
        win=21
        modeprofile=False
        seqs=readFasta("tests/ceba_bovin.fasta")
        modefasta=True
        labels='F'
        params=Params(weight,win)
        output = StringIO.StringIO()
        for name in seqs.keys():
            seqr = SeqRecord(seqs[name], id=name)
            seqr = pred_coil(seqr, params)
            printCoil(seqr,labels,io=output)
        bench_output = open("tests/ceba_bovin.Coils").read()
        self.assertEqual(output.getvalue(),bench_output)        

if __name__ == "__main__":
    import unittest
    unittest.main()
