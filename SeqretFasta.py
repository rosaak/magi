from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_protein

class SeqRet():
	'''
	Class to is a is same as seqret in emboss
	Requires Bio SeqIO
	'''

	def __init__(self,args):
		'''
		:param filename: Give only file not the file handle
		'''
		# Verify the input types
		if (isinstance(args,list)):
			self.filelist = args
		else:
			self.fh = args

	def getFastaRecord(self,recDict,seq_name):
		'''
		BioSeq fasta record dict -> single record
		'''
		return recDict.get_raw(str(seq_name)).decode()

	def getSingleFasta(self,seq_name):
		'''
		:param seq_name:
		:return: single fasta sequence
		'''
		return self.makeFastaDict().get_raw(str(seq_name)).decode()

	def lenRecord(self,recDict):
		'''
		:return: the length of Record
		'''
		return len(recDict)

	def makeFastaDict(self):
		'''
		:return: gives the index as hash if duplicate sequence then it gives error
		'''
		try :
			record_dict = SeqIO.index(self.fh, "fasta")
			return record_dict
		except ValueError as e:
			print (" Remove the duplicate sequences from file {}".format(e))
			exit(-1)
	"""		
	def makeFastaDict2(self):
		'''
		:return: gives the Sqlite index 
		'''
		
		idx_name = 'dbIndex.idx'
		try :
			r = SeqIO.index_db( idx_name, self.filelist, "fasta", generic_dna )
			return r
		except ValueError as e:
			print (" Remove the duplicate sequences from file {}".format(e))
			#exit(-1)
		#except :
                 #   pass
        """
                  
if __name__== '__main__':
	pass
	
