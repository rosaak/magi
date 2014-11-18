__author__ = 'Roshan Padmanabhan'
__version__= '0.3'

'''
This script takes up a config file and parse it.
Checks the nuc and pep files 
    - get the number of sequences in each files
    - checks whether any duplicate sequences are there 
        - if there then the program aborts
    - store these information in a log file
Runs the ProteinOrtho program
'''


import os,argparse
import shlex, subprocess
import fileinput
from SeqretFasta import SeqRet

class diskwalk(object):
    """Class for getting directory walking collections"""
    def __init__(self, path):
        self.path = path.rstrip('/')
    
    def dirExists(self):
        return os.path.exists(self.path)

    def enumeratePaths(self):
        """Returns the path to all the files in a directory as a list"""
        path_collection = []
        for dirpath, dirnames, filenames in os.walk(self.path):
            for file in filenames:
                fullpath = os.path.join(dirpath, file)
                try:
                    if os.path.exists(fullpath) :
                        path_collection.append(fullpath)
                except Exception as e:
                    raise e
                    exit(1)
        return path_collection

    def enumerateFiles(self):
        """Returns all the files in a directory as a list"""
        file_collection = []
        for dirpath, dirnames, filenames in os.walk(self.path):
            for file in filenames:
                file_collection.append(file)
        return file_collection

    def enumerateDir(self):
        """Returns all the directories in a directory as a list"""
        dir_collection = []
        for dirpath, dirnames, filenames in os.walk(self.path):
            for dir in dirnames:
                dir_collection.append(dir)
        return dir_collection

    def checkfile(self,fname):
        self.fname = fname
        if os.path.isfile(self.fname) is False:
            print ("\nError : File don't exists\n")
        else :
            return os.path.isfile(self.fname)

    def fullFilePath(self,fname):
            self.fname = fname
            return os.path.join(self.path,self.fname)

    def sepNucs(self):
            '''(list)-->list
            list out the full path of all the nuc files in the folder
            '''
            self.filelist=self.enumerateFiles()
            nuclist=[]
            [nuclist.append(self.fullFilePath(i)) for i in self.filelist if '.NUC' in i.upper()]
            return nuclist
            

    def sepPeps(self):
            '''(list)-->list
            list out the full path of all the pep files in the folder
            '''
            self.filelist=self.enumerateFiles()
            peplist=[]
            [peplist.append(self.fullFilePath(i)) for i in self.filelist if '.PEP' in i.upper()]
            return peplist
            
class mergeFiles(object):
    """docstring for mergeFiles"""
    def __init__(self, args):
        super(mergeFiles, self).__init__()
        self.fl = args
    
    def mergeflist(self):
        with open('/tmp/joinedNAFasta.fasta','w') as fout:
            try :
                for line in fileinput.input(self.fl):
                    fout.writelines(line)
            except:
                print("problem saving file.")
        return '/tmp/joinedNAFasta.fasta'




    #def fileExists(self):
    #   collect_truth=[]
    #   for files in self.enumeratePaths():
    #       retutn os.path.exists(files)
            
class config2Data(object):
    """This Class investigates contents of config file"""
    def __init__(self, dik):
        super(config2Data, self).__init__()
        self.v=['org_name','in_path','out_path','res_name']
        #print(dik.keys())
        #print(dik.values())
        #print(len(dik))
        if len(dik) < 3:
            print("Check config file\n")
            exit(-1)
        if len(dik) == 4 :
            self.lofg=dik[self.v[0]]
            self.inp=dik[self.v[1]]
            self.oup=dik[self.v[2]]
            self.rname=dik[self.v[3]]
        else:
            self.lofg=dik[self.v[0]]
            self.inp=dik[self.v[1]]
            self.oup=dik[self.v[2]]
    
    def check_seq_files(self):
        return self.lofg


    def check_number_of_files(self,l,p):
        '''(list,path)->bool
        
        >>>check_number_of_files(lofg,path)
        ...True
        >>>check_number_of_files(lofg,path)
        ...False
        
        Returns bool if the path contains files which is double 
        the number in list of genomes 
        
        use class method diskwalk.enumerateFiles New Empty File

        '''
        self.lofg=l
        self.inp =p
        d = diskwalk(self.inp)
        return len(d.enumerateFiles()) == len(self.lofg) *2

class checkConfigFiles(object):
    """This class checks the config files and paths
    method checkfile  : checks whether the file exits returns bool
    method conf2dict  : returns the dict from 
    """
    def __init__(self, fname,res_name):
        super(checkConfigFiles, self).__init__()
        self.fname = fname
        self.res_name = res_name
        if self.res_name:
            self.res_name = 'res_name'
            #print("The Result name file is : ",self.res_name)
        self.v=['org_name','in_path','out_path','res_name']


    def checkfile(self):
        self.fname
        if os.path.isfile(self.fname) is False:
            print ("\nError : File don't exists\n")
        else :
            return os.path.isfile(self.fname)
        
    def conf2dict(self):
        spl = [ lines.strip()  for lines in open(self.fname,'r').readlines()]
        pdict={}
        
        for l in spl:
            try:
                if l.startswith(self.v[0]):
                    names = (l.split(':')[1].split(','))
                    names = [i.strip() for i in names]
                    pdict[self.v[0]]=names
                elif l.startswith(self.v[1]):
                    path = (l.split(':')[1].strip())
                    pdict[self.v[1]]=path
                elif l.startswith(self.v[2]):
                    opath = (l.split(':')[1].strip())
                    pdict[self.v[2]]=opath
                #elif l.startswith(self.v[3]):
                #   r = (l.split(':')[1].strip())
                #   pdict[self.v[3]]=r
                else :
                    pass
            except KeyError as e :
                print ("Please check the config file again",e)
            if self.res_name  and l.startswith(self.v[3]):
                r = (l.split(':')[1].strip())
                pdict[self.v[3]]=r
        return pdict

class ProteinOrthoRun(object):
    """ Methods for running ProteinOrthoRun and """
    def __init__(self, genome_list, output_dir,results_file):
        super(ProteinOrthoRun, self).__init__()
        self.genome = genome_list # list containg genome file names
        #base= self.genome.split('/')[-1].split('.')[0]
        #self.nuc = base + ".nuc"
        #self.pep = base + ".pep"
        self.odir = output_dir.rstrip('/')
        self.rfile = results_file

    def P_cmdline(self):
        '''
        >>>ProteinOrthoRun.P_cmdline()

        '''
        #cmdline = ()
        #Proteinortho.pl blast=blastp E-value=1e-05, alg.-conn.=0.5, coverage=0.5, percent_identity=25, adaptive_similarity=0.95, inc_pairs=0, inc_singles=0, selfblast=0, unambiguous=1
        #proteinortho5.pl -p=blastp -conn=0.5 -identity=30 
        #cmd = "prodigal -i " + self.genome + " -o -c -m -q -g " + self.gcode + " -f sco  -a " + self.pep + " -d " + self.nuc 
        cmd = "proteinortho5.pl " + " -p=blastp " + " -conn=0.5 " + " -identity=30 " + ' '.join(self.genome) +" "+ " >"+self.odir + "/" + self.rfile
        c =  shlex.split(cmd)
        return c

    def P_run(self,args):
        '''
        
        Returns two files nuc and pep

        >>>ProteinOrthoRun.P_run()
        True
        >>>ProteinOrthoRun.P_run()
        False

        '''
        self.args = args
        try :
            s = subprocess.call(args, shell=False)
        except :
            raise
        return s

        
#######################         MAIN        ##########################

#######################         Command Line Argument Parsing       ##########################

if __name__ == '__main__':

    des="""
    This script runs ProteinOrtho which will be later used in AGIOs analysis.\n
    Requirements:\tProteinOrtho , blastall , mcl.\n
    """
    parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter , conflict_handler='resolve')
    parser.add_argument('-g', help='path to the genome files',nargs='+',  action='store',dest='genome_list',required=True)
    parser.add_argument('-o', help='out put directory path', action='store',dest='out_path',required=True )
    parser.add_argument('-r', help='results file name', action='store',dest='res_name')#,required=False)
    parser.add_argument('-x', help='Execute ProteinOrtho', action='store_true')
    args = parser.parse_args()

    glist = args.genome_list
    output_dir = args.out_path
    results_file = args.res_name
    pOrtho = args.x

    if  (pOrtho) and (not results_file):
        print("\n=====Error=====\nparameter -r is required if -x is on\n")
        parser.print_help()
        exit(1)
    elif  (not pOrtho) and (results_file):
        #print(pOrtho,results_file)
        pOrtho = True

    print("The parsed arguments are :-\n\n\
        the genome list : {}\n\
        the output directory : {}\n\
        the results file name : {}\n\
        the argument if protein ortho has to run: {}\n".format(glist,output_dir,results_file,pOrtho))


    # Create an object for checking the configfile

    if pOrtho:
        n = ProteinOrthoRun(glist,output_dir,results_file)
        print(n.P_cmdline())



