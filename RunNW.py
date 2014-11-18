#!/usr/bin/env python3

__author__ = 'Roshan Padmanabhan'
__version__ = '0.01'


import os , re, shlex, subprocess, argparse 
from SeqretFasta import SeqRet

class EmbossRun(object):
    """Class EmbossRun creates the instance of either needle or streacher
    The EmbossCmdline creates a 
    """
    def __init__(self, arg):
        super(EmbossRun, self).__init__()
        """
        seq is a BioSeq object
        """
        self.seq1 = arg[0]
        self.seq2 = arg[1]
        self.algo = arg[2] 
        self.emboss_out = 'emboss_out'

    def EmbossCline2(self):
        if self.algo == 'needle' :
            emboss_cmd = "needle " + " asequence=" + self.seq1 + " bsequence=" + self.seq2 + " gapopen=" + '10' + " gapextend=" + '0.5' + " outfile=" + self.emboss_out + " datafile=" + "EDNAFULL " + " aformat=" + "srspair " + "-auto"
            c =  shlex.split(emboss_cmd)
        else :
            emboss_cmd = "stretcher " + " asequence=" + self.seq1 + " bsequence=" + self.seq2 + " gapopen=" + '16' + " gapextend=" + '4' + " outfile=" + self.emboss_out + " datafile=" + "EDNAFULL " + " aformat=" + "srspair " + "-auto"
            c =  shlex.split(emboss_cmd)
        return c
    
    def runCline(self,args):
        self.args =args
        try :
            s = subprocess.call(args, shell=False)
        except :
            raise
        return s

    def GetSimilarity(self):
        '''
        Returns the similarity the needle or streacher file  
        '''
        N=open(self.emboss_out,'r').readlines()
        sim = re.compile('# Similarity:')
        for m in N:
            try:
                if (sim.match(m)):
                    #print(m.split('\n')[0].split("(")[1].split("%")[0])
                    return m.split('\n')[0].split("(")[1].split("%")[0]
            except:
                return 'NA'
                continue

    def GetIdentity(self):
        '''
        Returns the identity the needle or streacher file  
        '''
        N=open(self.emboss_out,'r').readlines()
        ide = re.compile('# Identity:')
        for m in N:
            try:
                if (ide.match(m)):
                    #print (m.split('\n')[0].split("(")[1].split("%")[0])
                    return m.split('\n')[0].split("(")[1].split("%")[0]
            except:
                return 'NA'
                continue

if __name__ == '__main__':
    fl = ['.agiosseq1.fasta','.agiosseq2.fasta','needle']
    #fl = ['.agiosseq1.fasta','.agiosseq2.fasta','stretcher']
    x = EmbossRun(fl)
    print(x.seq1,x.seq2,x.algo,x.emboss_out)
    print(x.EmbossCline2())
    ncl =x.EmbossCline2()
    x.runCline(ncl)
    print( x.GetSimilarity()) 
    print( x.GetIdentity()) 

