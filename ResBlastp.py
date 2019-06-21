import os,csv,sys
import FileHandle as fh

#this part of the module is used to run blast between genomes
#and databases of antibiotic resistance genes
class ResBlastp:
	#attributes of class Blast
	#including software install directory, fasta files directory
	#blast files directory, name list of cds fasta files and blast results
	install_dir = ""
	in_dir = ""
	out_dir = ""
	blast_dir = ""
	fasta_files = []
	blast_files = []

	def __init__(self,install_directory,in_directory,out_directory):
		#get the install directory of the software
		self.install_dir = fh.cwd_get(install_directory)
		#get the cds files directory
		self.in_dir = fh.inpd_get(in_directory,file_categories="faa files")
		self.out_dir = fh.inpd_get(out_directory,file_categories="ouput files")
		#get all the cds file names form cds_directory
		self.fasta_files = fh.filename_get(self.in_dir,".faa")
		if len(self.fasta_files) == 0:
			print("No fasta file was found, stop analyzing...")
			sys.exit(0)
		#get the install directory of blast+
		#according to the setting modified by the user
		self.blast_dir = fh.setting_reader(self.install_dir,"blast+")
		#get all the blast results exist
		self.blast_files = fh.filename_get(self.out_dir,"res_prot.xls")

	#delete all the blank and tab in the fasta files
	#to avoid find wrong contigs in the later module
	def blankdel(self):
		for each_name in self.fasta_files:
			content = []
			each = self.in_dir+each_name+".faa"
			with open(each,"r") as f:
				content = f.read()
			content = content.replace(" ","_")
			content = content.replace("	","|")
			content = content.replace(",","_")

			with open(each,"w") as f:
				f.write(content)

	def makedb(self):
		print("\nstart making blast database...")
		#get the databases from the install directory
		db_name = self.install_dir+"databases/res_prot"
		#create blast databases into the blast+ directory
		makeblastdb = self.blast_dir+"makeblastdb"
		os.system(makeblastdb+" -in "+db_name+" -out res_prot -dbtype prot")
		print("\nfinish making blast database...\n")

	def blastp(self):

		def os_blast(in_name,query_name,output_name):
			print(in_name+"	blast start!!!")
			#add the directory of blastn
			blastp = self.blast_dir+"blastp"
			#run blastn in the terminal with parameters according to blast+ manual
			os.system(blastp+" -query "+query_name+" -db res_prot -out "+output_name+\
				" -outfmt \"6 std slen\" -evalue 1e-20 -num_alignments 6")
			print(in_name+"	blast finished!!!")

		#if blast results exist, you can choose not to re-blast them,
		#because it may waste a lot of time.
		#but if you change the input files, please update your blast results.
		if len(self.blast_files) == len(self.fasta_files):
			#you can choose to update them one by one, all or not to
			update_choose = input("update all blast results (O(one by one)/A(all)/N(no))?:").upper()
			for each_name in self.fasta_files:
				#get the fasta file name of the genome
				each = self.in_dir+each_name+".faa"
				#get the output file name of blast
				out_name = self.out_dir+each_name+"VSres_prot.xls"
				#working when choose to blast one by one
				if update_choose == "O":
					#choose to blast or not to
					choose = input("update blast result for "+each_name+" (Y/N)?:").upper()
					if choose == "Y":
						fh.file_del(out_name)
						os_blast(each_name,each,out_name)
					else:
						print(each_name+"	blast skipped...")
				#working when choose to re-blast all genomes
				elif update_choose == "A":
					fh.file_del(out_name)
					os_blast(each_name,each,out_name)
				#skipping when choose not to re-blast any of them
				else:
					print(each_name+"	blast skipped...")
		#if the blast results are not equivalent to the input files,
		#blast will start from the beginning to the end
		else:
			for each_name in self.fasta_files:
				each = self.in_dir+each_name+".faa"
				out_name = self.out_dir+each_name+"VSres_prot.xls"
				fh.file_del(out_name)
				os_blast(each_name,each,out_name)

#this part of the module is used to parse the result of blast
#annotation files are written into "_ar.csv"
class ParseResblastResult:
	#attributes of class ParseBlastResult
	#including software install directory, blast result files directory
	#name list of blast results and annotation files
	#identity and coverage threshold while parsing blast results
	install_dir = ""
	in_dir = ""
	out_dir = ""
	arg_dir = ""
	blast_files = []
	annotation_files = []
	identity_p = 98.0
	query_coverage_p = 0.98

	def __init__(self,install_directory,in_directory,out_directory):
		#get the install directory of the software
		self.install_dir = fh.cwd_get(install_directory)
		#get the blast result files directory
		self.in_dir = fh.inpd_get(in_directory,file_categories="blast files")
		self.out_dir = fh.inpd_get(out_directory,file_categories="ouput files")
		#create the new directory to save files of re-annotated ar gene sequences
		self.arg_dir = self.out_dir+"arg/"
		fh.dir_add(self.arg_dir)
		#add names of blast result files into a list
		self.blast_files = fh.filename_get(self.out_dir,"res_prot.xls",showext=True)
		#get the threshold of identity
		self.identity_p = float(fh.setting_reader(self.install_dir,"identity_p"))
		#get the threshold of query coverage
		self.query_coverage_p = float(fh.setting_reader(self.install_dir,"query_coverage_p"))

	#parse all the blast result files in the list "blast_files"
	#return annotation files for each genome in the same directory
	def parse(self):

		#annotation method for each blast file
		def annotation(blast_file_name):
			print(blast_file_name+"	parse start...")
			#get the genome name from blast file name
			each_name = blast_file_name.split("VS")[0]
			#create the output annotation file name
			ar_csv = self.out_dir+each_name+"_ar.csv"
			#add the annotation file name into a list
			self.annotation_files.append(each_name+"_ar.csv")
			fh.file_del(ar_csv)
			#read the blast result and store into a list
			blast_content = fh.file_reader(self.out_dir+blast_file_name)
			#create the annotation file
			f_ar = open(ar_csv,"a")
			#write the title of annotation file
			#the origin, start and end of predicted ar genes are added compared to blastp
			f_ar.write("#gene_num"+","+"gene_name"+","+"identity"+","+"e_value"\
				+","+"ARO_num"+","+"accession_num"+","+"ar_gene_allele"+","+"drug_class"\
				+","+"ar_mechanism"+"\n")	
			gene_identity = {}
			gene_perc = {}
			write_content = {}

			#get the information of the gene according to ARO
			#return the gene family, drug class and ar mechanism of the gene
			def nucl_arinfo_get(ar_info_list):
				str_ar_gene = ""
				str_drug_class = ""
				str_ar_mechan = ""
				ar_gene_allele_r = []
				drug_class_r = []
				ar_mechanism_r = []

				for row in annotation_info:
					#find the information according to ARO
					#append the information into lists
					if ar_info_list[0] == row[0]:
						str_ar_gene = row[-1]
						str_drug_class = row[1]
						str_ar_mechan = "/"

				return str_ar_gene,str_drug_class,str_ar_mechan

			#deal with each row in the blast result
			for line in blast_content:
				#separate each item into a list
				items = line.split("	")
				#separate subject sequence information into a list
				ar_info = items[1].split("_")
				#get the detailed infomation from the database
				#including antibiotic gene family, drug classes and mechanisms
				ar_gene_allele,drug_class,ar_mechanism = nucl_arinfo_get(ar_info)

				#only be annotated when higher than the threshold of coverage
				perc = float(items[3])/float(items[12])
				if float(items[2]) >= self.identity_p and perc > self.query_coverage_p:
					#if the first gene in a new scaffold
					if items[0] not in gene_identity.keys():
						#append the gene identity into a dictionary
						gene_identity[items[0]] = float(items[2])
						gene_perc[items[0]] = perc
						#append the content that will be written into a dictionary
						write_content[items[0]] = items[0]+","+items[1]+","+items[2]+\
						","+items[10]+",/,/,"+ar_gene_allele+","+drug_class+","+ar_mechanism+"\n"
					else:
						if float(items[2]) > gene_identity[items[0]]:
							gene_identity[items[0]] = float(items[2])
							gene_perc[items[0]] = perc
							write_content[items[0]] = items[0]+","+items[1]+","+items[2]+\
							","+items[10]+",/,/,"+ar_gene_allele+","+drug_class+","+ar_mechanism+"\n"
						elif float(items[2]) == gene_identity[items[0]] and perc > gene_perc[items[0]]:
							gene_perc[items[0]] = perc
							write_content[items[0]] = items[0]+","+items[1]+","+items[2]+\
							","+items[10]+",/,/,"+ar_gene_allele+","+drug_class+","+ar_mechanism+"\n"


			#write the annotation content into the annotation file
			for key in write_content.keys():
				f_ar.write(write_content[key])
			print(blast_file_name+"	parse finished!!!")
			f_ar.close()

		#read the databases informations and store into a list
		annotation_info = fh.csv_reader(self.install_dir+"databases/res_parse")
		for file in self.blast_files:
			annotation(file)

	def ar_write(self):

		def each_write(csvfilename,oldextension,newextension):
			ar_genes = {}
			start_end = {}
			#get the old cds file name of genome
			old_name = self.in_dir+genome_id+oldextension
			#create the new cds file name of genome
			new_name = self.arg_dir+genome_id+"_ar"+newextension
			f_new = open(new_name,"w")
			#read the old cds file content and store in a list
			old_content_r = fh.file_reader(old_name)
			old_content = []
			#get the index of >scaffold in the file
			scaf_index = [old_content_r.index(line) for line in old_content_r if ">" in line]
			scaf_index.append(len(old_content_r))
			#to locate the gene in a whole scaffold
			#gather a scaffold which may be separated by "\n"
			for i in range(len(scaf_index)-1):
				scaf_each = ""
				#title of the scaffold
				old_content.append(old_content_r[scaf_index[i]])
				#all the sequences of the scaffold
				for x in range(scaf_index[i]+1,scaf_index[i+1]):
					scaf_each += old_content_r[x].strip("\n")
				old_content.append(scaf_each)

			new_content = []
			#read related csv annotation file and store arg information in a dictionary
			plots = fh.csv_reader(self.out_dir+csvfilename)
			for row in plots:
				ar_genes[row[0]] = " ["+row[1]+"|identity:"+row[2]+"|"+row[4]+"|"+row[5]+"]"
			for gene in ar_genes.keys():
				for i in range(len(old_content)):
					if ">" in old_content[i]:
						#add new ar gene name and seq into new_content list
						if gene in old_content[i]:
							new_content.append(">"+gene+ar_genes[gene]+"\n")
							new_content.append(old_content[i+1]+"\n")
			for line in new_content:
				f_new.write(line)
			f_new.close()

		#for each genome, invoke fuction each_write
		fna_count = 0
		for filename in os.listdir(self.in_dir):
			if filename.endswith(".fna"):
				fna_count += 1
		fna_extract = False
		if fna_count == len(self.annotation_files):
			fna_extract = True

		#for each genome, invoke fuction each_write
		for each in self.annotation_files:
			genome_id = each.split("_ar.csv")[0]
			each_write(each,".faa",".faa")
			if fna_extract:
				each_write(each,".fna",".fna")
			print("finish writing ar sequences of "+genome_id)


def main(install_directory,in_directory,out_directory):
	resblast = ResBlastp(install_directory,in_directory,out_directory)
	resblast.blankdel()
	resblast.makedb()
	resblast.blastp()
	pbr = ParseResblastResult(resblast.install_dir,resblast.in_dir,resblast.out_dir)
	pbr.parse()
	pbr.ar_write()
	return resblast.out_dir

if __name__ == '__main__':
	main("","","")