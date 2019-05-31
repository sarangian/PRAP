from Bio import SeqIO
import FileHandle as fh
import csv,os,sys

#this part is used to extract pretein coding sequence from
#the genebank annotation files, protein sequences in fasta format
#and nucleotide sequences in fasta fomat will be generated

class CdsExtract:
	#attrbutes of class CdsExtract
	#including the genebank file directory and name list
	gb_file_dir = ""
	gb_genome_id = []

	def __init__(self,file_directory):
		#get the directory of genebank files
		self.gb_file_dir = fh.inpd_get(file_directory,file_categories="genebank files")
		#change work directory to genebank files
		os.chdir(self.gb_file_dir)
		#get all the gb files into a list
		self.gb_genome_id = fh.filename_get(self.gb_file_dir,".gb")
		if len(self.gb_genome_id) == 0:
			print("No genebank file was found, stop analyzing...")
			sys.exit(0)

	def cds_extract(self):
		for each in self.gb_genome_id:
			print("Extracting CDS from "+each)
			filename = each+".gb"
			#form the fasta format filename of the handling genome
			cds_prot_file = each+".faa"
			cds_nucl_file = each+".fna"
			fh.file_del(cds_prot_file)
			fh.file_del(cds_nucl_file)
			f_cds_prot = open(cds_prot_file,"a")
			f_cds_nucl = open(cds_nucl_file,"a")
			#extract cds from genebank files
			gene_count = 0
			records = list(SeqIO.parse(filename, "genbank"))
			for record in records:
				for feature in record.features:
					#only write the sequences with type "CDS" and with translation sequences
					if feature.type == "CDS" and "translation" in feature.qualifiers:
						#give each gene a new name in order like gene_1, gene_2,...geng_m
						gene_count += 1
						f_cds_nucl.write(">"+"gene_"+str(gene_count)+"|"+each+"\n")
						#the product of the gene is also recorded after the gene name like [transportase]
						f_cds_nucl.write(str(record.seq)[int(feature.location.start):int(feature.location.end)]+"\n")
						#write the nucleotide sequence
						f_cds_prot.write(">"+"gene_"+str(gene_count)+"|"+each+"\n")
						#write the same gene name of protein sequence with nucleotide sequence
						f_cds_prot.write(feature.qualifiers["translation"][0]+"\n")
						#write the protein sequence
			f_cds_nucl.close()
			f_cds_prot.close()

def main(file_directory):
	cds = CdsExtract(file_directory)
	cds.cds_extract()
	return cds.gb_file_dir

if __name__ == '__main__':
	#run main() when use this module separately
	main("")