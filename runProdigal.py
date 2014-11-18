__author__ = 'Roshan'
__version__ = '0.1'


import shlex
import subprocess
import argparse
import os 


class Run_Prodigal(object):
	""" Methods for running Prodigal and gives an edited version of fasta files """
	def __init__(self, genomefilename):
		super(Run_Prodigal, self).__init__()
		self.genome = genomefilename # single genome file
		self.gcode = '11'
		base= self.genome.split('/')[-1].split('.')[0]
		self.nuc = base + ".nuc"
		self.pep = base + ".pep"
		self.sco = base + ".sco"

	def P_cmdline(self, gcode=None):
		'''
		Returns str of  prodigal commandline
		'''
		if gcode is None:
			self.gcode == str(gcode)
		else :
			self.gcode == str(gcode)

		cmd = "prodigal -i " + self.genome + " -c -m -q -g " + self.gcode + " -o sco  -a " + self.pep + " -d " + self.nuc 
		c =  shlex.split(cmd)
		return c


	def P_run(self,args):
		'''
		
		Returns True if prodigal runs

		>>>Run_Prodigal.P_run()
		True
		>>>Run_Prodigal.P_run()
		False

		'''
		self.args = args
		try :
			s = subprocess.call(args, shell=False)
		except :
			raise
		return s


	def shorten_fasta_header(self,fn):
		'''
		creates another edited version of input fasta file
		'''
		with open(fn,'r') as fh:
			fr = fh.readlines()
			nfn = fn + ".modified"
			fw = open(nfn,'w+')
			for i in fr:
				if i.startswith(">"):
					fw.write(i.split(" ")[0]+"\n")
				else:
					fw.write(i)
			fw.close()
		return nfn

	def remove_files(self):
		'''
		just remove the files from the disk
		'''
		os.remove(self.nuc)
		os.remove(self.pep)
		#os.remove(self.sco)
	def rename_files(self):
		'''
		just remove the files from the disk
		'''
		nucfile = self.nuc + ".modified"
		pepfile = self.pep + ".modified"
		os.rename(nucfile,self.nuc)
		os.rename(pepfile,self.pep)
		



if __name__ == '__main__':
	
	######################			Command Line Argument Parsing  		##########################
	des="""
	\n
	This is a suplimentary script which runs prodigal and parse the file and give correct filenames\n
	Requirements: prodigal\n
	"""

	parser = argparse.ArgumentParser(description=des,formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-g', help='path to the genome files',nargs='+', action='store',dest='genome_list',required=True)
	args = parser.parse_args()
	x = args.genome_list

	for i in x:
		#checkfile(i)
		n = Run_Prodigal(i)
		cln = n.P_cmdline()
		print("running prodigal for {}".format(n.genome))
		n.P_run(cln)
		n.shorten_fasta_header(n.nuc)
		n.shorten_fasta_header(n.pep)
		n.remove_files()
		n.rename_files()
		print("chek result files for {}  {} ".format(n.nuc, n.pep))
	

'''

TODO : Checking Files Before Running Should Be Included 

def checkfile(fname):
	if os.path.isfile(fname) is False:
		print ("\nError : File don't exists\n")
	else :
		return os.path.isfile(fname)

'''