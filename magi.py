__author__ = '''Roshan Padmanabhan'''
__version__ = '0.1'


'''
Main MAGI program
Absolute Requirement : pandas , numpy , biopython all compatabile with python3 
Other Programs Required to run : Emboss , ProteinOrtho, Blastall for running ProteinOrtho , Prodigal
Input is a directory containing all pep and nuc files
pep and nuc files can be generated using another script called runProdigal.py 


'''


#########################################################
import re
import sys
import os
import argparse
import shlex
import subprocess
from itertools import combinations as comb
import pandas as pd
import numpy as np
from SeqretFasta import SeqRet 
import SeqretFasta
from preProcess import diskwalk , mergeFiles
import RunNW
#########################################################

def nC2(list_containing_headers):
    '''list->list
    list_containg_headers has the header names
    returns the list of nC2 combiantions in tuple
    '''
    Comb = comb(list_containing_headers,2)
    List_of_Tuples = []
    for i in Comb:
        List_of_Tuples.append(i)
    return List_of_Tuples

def Mmax(l):
    '''
    return the max value of the two
    '''
    if len(l) == 2:
        return max(l)
    else :
        return "Problem"

def MakeNewDF(df1):
    '''(df) -> (df)
    cut the first three columns and
    returns the rest of columns as DF
    '''
    lookup1= ['#species','proteins','alg.-conn.']
    df2 = df1.drop(lookup1, 1)
    df2 = df2.drop(df2.index[:1])
    return df2

def getCore(df):
    '''(df) -> (df)
    returns the core set in all rows
    input is the dataframe returned from MakeNewDF
    it reomoves all rows with '*'
    '''
    mask = df.applymap(lambda x: x in ['*'])
    newDF = df[-mask.any(axis=1)]
    return newDF

def getHeadNames(pdDF):
    '''df->list
    Returns the headers of pandas data fram as a list
    >>>getHeadNames(DF1)
    ['sp1', 'sp2', 'sp3']
    '''
    return [ hNames for hNames in  pdDF.columns]

def ConcatHeaders(lst):
    '''
    list->list
    returns the concatenated list of Species
    First it does a n combianation 2 of the names in list
    and then concatenate the combinations using '_vs_'
    then returns the list
    >>>ConcatHeaders(['a','b','c','d'])
    ['a_vs_b', 'a_vs_c', 'a_vs_d', 'b_vs_c', 'b_vs_d', 'c_vs_d']
    >>>ConcatHeaders(['1','2','3','4'])
    ['1_vs_2', '1_vs_3', '1_vs_4', '2_vs_3', '2_vs_4', '3_vs_4']
    >>>ConcatHeaders([1,2,3,4])
    ['1_vs_2', '1_vs_3', '1_vs_4', '2_vs_3', '2_vs_4', '3_vs_4']
    '''
    return [str(x[0])+"_vs_"+str(x[1]) for x in nC2(lst)]

def ConcatIndex(lst):
    '''
    list->list
    returns the concatenated list of seqids using '_vs_'
    then returns the list
    >>>ConcatIndex(['seqid1', 'seqid2'])
    ['seqid1_vs_seqid2']
    >>>ConcatIndex([1,2])
    ['1_vs_2']
    '''
    return [str(lst[0])+"_vs_"+str(lst[1])]

def ConcatIndexString(lst):
    '''
    list->str
    returns the concatenated str of seqids using '_vs_'
    then returns the list
    >>>ConcatIndex(['seqid1', 'seqid2'])
    'seqid1_vs_seqid2'
    >>>ConcatIndex([1,2])
    '1_vs_2'
    '''
    return str(lst[0])+"_vs_"+str(lst[1])

def running_emboss(emboss_parameters):
    '''(list) -> (str)
    returns the similarity values for the parameters provided
    emboss_parameters is a list comtaining seq1filename , seq2filename and needle / streacher
    '''
    erun = RunNW.EmbossRun(emboss_parameters)
    ncl = erun.EmbossCline2()
    erun.runCline(ncl)
    return erun.GetSimilarity()

def create_pSeries(data,name):
    '''
    returns a pandas series
    data should be dict and name string
    '''
    name =  pd.Series(data,index=None,name=name,dtype=float)
    return name

def create_pDataFrame(list_of_series,list_of_index_names):
    '''list,list->DataFrame
    returns a pandas dataframe
    '''
    return pd.DataFrame(list_of_series,index=list_of_index_names).T

def getMean(pdSeries):
	'''pandas series --> number
	'''
	if getLen(pdSeries) >= 1:
		return round(pdSeries.mean(),2)

def getMin(pdSeries):
	'''pandas series --> number
	'''
	return pdSeries.min()

def getMax(pdSeries):
	'''pandas series --> number
	'''
	return pdSeries.max()

def getLen(pdSeries):
	'''pandas series --> number
	'''
	return pdSeries.count()

def garbage(anything):
	del anything

def dots():
	print(".",end='')

def get_agios_for_two_dict(dfcolumn1,dfcolumn2):
    templistof_org1_org2 = {}
    for i,en in enumerate(range(1,dfcolumn1.count())):
        seq1outfile = open("/tmp/.agiosseq1.fasta",'w')
        seq2outfile = open("/tmp/.agiosseq2.fasta",'w')
        if '*' not in [ dfcolumn1[en] , dfcolumn2[en] ]:
            seq1outfile.write(s.getSingleFasta( dfcolumn1[en] ))
            seq2outfile.write(s.getSingleFasta( dfcolumn2[en] ))
            seq1outfile.close()
            seq2outfile.close()
            if not (len(SeqretFasta.SeqIO.read(seq1outfile.name,'fasta'))) > 1500 :
                templistof_org1_org2[(ConcatIndexString([dfcolumn1[en],dfcolumn2[en]]))] = running_emboss([seq1outfile.name,seq2outfile.name,'needle'])        
            else :
                templistof_org1_org2[(ConcatIndexString([dfcolumn1[en],dfcolumn2[en]]))] = running_emboss([seq1outfile.name,seq2outfile.name,'needle'])        
    return templistof_org1_org2

#########################################################

if __name__ == '__main__':

    #######################         Command Line Argument Parsing       ##########################
    des="""
    The main MAGI script.\n
    Inputs are ProteinOrtho result file, path to the cds files and output directory path \n
    """
    parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter )
    parser.add_argument('-g', help='path to input directory', action='store',dest='genome_files',required=True)
    parser.add_argument('-p', help='path to ProteinOrtho result file', action='store',dest='proteinOrtho_Result_file',type=argparse.FileType('r'),required=True) 
    parser.add_argument('-o', help='Out put directory path', action='store',dest='out_path',required=True )
    parser.add_argument('-c', help='calculate agios for core genes (0 for no cores 1 for yes)',type=int,choices=range(0,2),dest='core',required=False)

    args = parser.parse_args()

    glist = args.genome_files
    output_dir = args.out_path
    po_res_file = args.proteinOrtho_Result_file.name
    core_agios = args.core

    # ==========================TO DO==============================
    # check if the directory is correct
    # ==========================TO DO==============================
    print("The input nucleotide files are processing now")
    d = diskwalk(glist)
    try:
        allnuc = mergeFiles(d.sepNucs())
        nucsmerged = allnuc.mergeflist()
        #print("the path to the meged file ",nucsmerged)
    except Exception as e:
        raise e
    o = diskwalk(output_dir)
    outfile = o.fullFilePath('AGIOSresults.csv')


    ############################### START PROCESSING FILES IN THE INPUT FOLDER #########################################

    # Now start the module SeqretFasta.SeqRet and get the sequence from the file nucsmeged 
    s=SeqRet(nucsmerged)
    recd=s.makeFastaDict()

    ################################# START PROCESSING PROTEINORTHO RESULT FILE #######################################    
    pOrthoresfn = po_res_file
    # The initial dataframe
    pOResults = pd.read_csv(pOrthoresfn,sep='\t')

    # corrected working DF with or without core genes
    if core_agios :
    	pdres = getCore(MakeNewDF(pOResults))
    else:
    	pdres = MakeNewDF(pOResults)

    garbage(pOrthoresfn)

    Headers = getHeadNames(pdres)
    #print("The headers",Headers)
    header_groups = nC2(Headers)
    #print("The header groups",ConcatHeaders(Headers))

    # final df to hold the results
    finalDF = pd.DataFrame(columns=ConcatHeaders(Headers),index=None)

    # looping tru all the header combinations
    for OrgTup in header_groups[:2]:
        org1=OrgTup[0] # org1 is the column header
        org2=OrgTup[1]
        colOrg1 = pdres[org1] # colOrg1 is the contents of the column org1
        colOrg2 = pdres[org2]
        print("============== For the pair {} == {} ============".format(org1,org2 ))
        dots()
        pSname = ConcatIndexString(OrgTup)
        #print("Calculating AGIOS for the pair {} and {}".format(org1,org2 ))
        #make a series for each header_group
        eachSeries = create_pSeries(get_agios_for_two_dict(colOrg1,colOrg2),pSname)
        # put it in the final df
        finalDF[pSname] = [getMean(eachSeries),str(getLen(eachSeries))]
    print('')
# print the results 
FDF = pd.DataFrame(finalDF.T)
FDF.columns = ['agios_value','number_of_orthologs']
FDF.to_csv(outfile)
print("Please check the file {} for results".format(outfile))
print("")
