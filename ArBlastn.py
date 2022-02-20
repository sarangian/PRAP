#!/usr/bin/python3
import os,csv,sys
import FileHandle as fh

#this part of the module is used to run blast between genomes
#and databases of antibiotic resistance genes
class Blastn:
	#attributes of class Blast
	#including software install directory, fasta files directory
	#blast files directory, name list of cds fasta files and blast results
	install_dir = ""
	in_dir = ""
	blast_dir = ""
	out_dir = ""
	fasta_files = []
	blast_files = []
	identity_n = 90.0

	def __init__(self,install_directory,in_directory,out_directory):
		#get the install directory of the software
		self.install_dir = fh.cwd_get(install_directory)
		#get the fasta files directory
		self.in_dir = fh.inpd_get(in_directory,file_categories="fasta files")
		self.out_dir = fh.inpd_get(out_directory,file_categories="ouput files")
		#get all the fasta file names form cds_directory
		self.fasta_files = fh.filename_get(self.in_dir,".fna")
		if len(self.fasta_files) == 0:
			print("No fasta file was found, stop analyzing...")
			sys.exit(0)
		#get the install directory of blast+
		#according to the setting modified by the user
		self.blast_dir = fh.setting_reader(self.install_dir,"blast+")
		#get all the blast results exist
		self.blast_files = fh.filename_get(self.out_dir,"ar_nucl.xls",showext=True)
		#get the threshold of identity
		self.identity_n = float(fh.setting_reader(self.install_dir,"identity_n"))

	#delete all the blank and tab in the fasta files
	#to avoid find wrong contigs in the later module
	def blankdel(self):
		for each_name in self.fasta_files:
			content = []
			each = self.in_dir+each_name+".fna"
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
		db_name = self.install_dir+"databases/ar_nucl"
		#create blast databases into the blast+ directory
		makeblastdb = self.blast_dir+"makeblastdb"
		os.system(makeblastdb+" -in "+db_name+" -out ar_nucl -dbtype nucl")
		print("\nfinish making blast database...\n")

	def blastn(self):

		def os_blast(in_name,query_name,output_name):
			print(in_name+"	blast start!!!")
			#add the directory of blastn
			blastn = self.blast_dir+"blastn"
			#run blastn in the terminal with parameters according to blast+ manual
			os.system(blastn+" -query "+query_name+" -db ar_nucl -out "+output_name+\
				" -outfmt \"6 std slen\" -evalue 1e-20 -perc_identity "+str(self.identity_n))
			print(in_name+"	blast finished!!!")

		#if blast results exist, you can choose not to re-blast them,
		#because it may waste a lot of time.
		#but if you change the input files, please update your blast results.
		if len(self.blast_files) == len(self.fasta_files):
			#you can choose to update them one by one, all or not to
			update_choose = input("update all blast results (O(one by one)/A(all)/N(no))?:").upper()
			for each_name in self.fasta_files:
				#get the fasta file name of the genome
				each = self.in_dir+each_name+".fna"
				#get the output file name of blast
				out_name = self.out_dir+each_name+"VSar_nucl.xls"
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
				each = self.in_dir+each_name+".fna"
				out_name = self.out_dir+each_name+"VSar_nucl.xls"
				fh.file_del(out_name)
				os_blast(each_name,each,out_name)

#this part of the module is used to parse the result of blast
#annotation files are written into "_ar.csv"
class ParseBlastnResult:
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
	query_coverage_n = 0.80

	def __init__(self,install_directory,in_directory,out_directory):
		#get the install directory of the software
		self.install_dir = fh.cwd_get(install_directory)
		#get the blast result files directory
		self.in_dir = fh.inpd_get(in_directory,file_categories="fna files")
		self.out_dir = fh.inpd_get(out_directory,file_categories="ouput files")
		#create the new directory to save files of re-annotated ar gene sequences
		self.arg_dir = self.out_dir+"arg/"
		fh.dir_add(self.arg_dir)
		#add names of blast result files into a list
		self.blast_files = fh.filename_get(self.out_dir,"ar_nucl.xls",showext=True)
		self.annotation_files = fh.filename_get(self.out_dir,"_ar.csv",showext=True)
		#get the threshold of query coverage
		self.query_coverage_n = float(fh.setting_reader(self.install_dir,"query_coverage_n"))

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
				+","+"ar_mechanism"+","+"origin"+","+"start"+","+"end"+"\n")	
			scaf = {}
			gene_identity = {}
			gene_perc = {}
			write_content = {}
			start_end = {}
			gene_count = 0

			#to identify whether the gene have common regions with another one
			#return a bool value
			def in_query(point1,point2,start,end):
				if int(point2) <= int(start) or int(point1) >= int(end):
					return False
				else:
					return True

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
					if ar_info_list[4] == row[0]:
						ar_gene_allele_o = row[-1].split(";")
						ar_gene_allele_r += [each for each in ar_gene_allele_o \
						if each not in ar_gene_allele_r]
						drug_class_o = row[4].split(";")
						drug_class_r += [each for each in drug_class_o \
						if each not in drug_class_r]
						ar_mechanism_o = row[5].split(";")
						ar_mechanism_r += [each for each in ar_mechanism_o \
						if each not in ar_mechanism_r]

				#write the items in the lists into strings
				for i in range(len(ar_gene_allele_r)):
					if i < len(ar_gene_allele_r)-1:
						str_ar_gene += ar_gene_allele_r[i]+";"
					else:
						str_ar_gene += ar_gene_allele_r[i]
				for i in range(len(drug_class_r)):
					if i < len(drug_class_r)-1:
						str_drug_class += drug_class_r[i]+";"
					else:
						str_drug_class += drug_class_r[i]
				for i in range(len(ar_mechanism_r)):
					if i < len(ar_mechanism_r)-1:
						str_ar_mechan += ar_mechanism_r[i]+";"
					else:
						str_ar_mechan += ar_mechanism_r[i]

				return str_ar_gene,str_drug_class,str_ar_mechan

			#deal with each row in the blast result
			for line in blast_content:
				#separate each item into a list
				items = line.split("	")
				#separate subject sequence information into a list
				ar_info = items[1].split("|")
				#get the detailed infomation from the database
				#including antibiotic gene family, drug classes and mechanisms
				ar_gene_allele,drug_class,ar_mechanism = nucl_arinfo_get(ar_info)

				#get the start and end of the query
				if int(items[6]) < int(items[7]):
					start = items[6]
					end = items[7]
				elif int(items[7]) < int(items[6]):
					start = items[7]
					end = items[6]

				#only be annotated when higher than the threshold of coverage
				#no gapopen is permitted to avoid frameshift mutation
				perc = float(items[3])/float(items[12])
				if perc >= self.query_coverage_n and int(items[5]) <= 0:
					#if the first gene in a new scaffold
					if items[0] not in scaf.keys():
						gene_count += 1
						#give a new name to the gene
						gene_num = each_name+"_AMRgene_"+str(gene_count)+"_from_"+items[0]
						scaf[items[0]] = [gene_num]
						#append the gene identity into a dictionary
						gene_identity[gene_num] = float(items[2])
						gene_perc[gene_num] = perc
						#append the start and the end of the gene into a dictionary
						start_end[gene_num] = [start,end]
						write_content[gene_num] = gene_num+","+ar_info[-1]+","+items[2]+\
						","+items[10]+","+ar_info[4]+","+ar_info[1]+","+ar_gene_allele+\
						","+drug_class+","+ar_mechanism+","+items[0]+","+start+","+end+"\n"
					#if not the first gene in the scaffold
					else:
						contain = False
						for each_gene in scaf[items[0]]:
							#identify whether the gene is the same with the former one
							contain = in_query(start,end,start_end[each_gene][0],start_end[each_gene][1])
							if contain:
								gene_num = each_gene
								break
						#if they have common regions
						if contain:
							#if the identity is higher than the former one
							#replace the former one with the latter one
							if float(items[2]) > gene_identity[gene_num]:
								gene_identity[gene_num] = float(items[2])
								gene_perc[gene_num] = perc
								start_end[gene_num] = [start,end]
								write_content[gene_num] = gene_num+","+ar_info[-1]+","+items[2]+\
								","+items[10]+","+ar_info[4]+","+ar_info[1]+","+ar_gene_allele+\
								","+drug_class+","+ar_mechanism+","+items[0]+","+start+","+end+"\n"
							elif float(items[2]) == gene_identity[gene_num] and perc > gene_perc[gene_num]:
								gene_identity[gene_num] = float(items[2])
								gene_perc[gene_num] = perc
								start_end[gene_num] = [start,end]
								write_content[gene_num] = gene_num+","+ar_info[-1]+","+items[2]+\
								","+items[10]+","+ar_info[4]+","+ar_info[1]+","+ar_gene_allele+\
								","+drug_class+","+ar_mechanism+","+items[0]+","+start+","+end+"\n"
						#if they don't have common regions
						elif not contain:
							#create a new gene name for a new ar gene
							gene_count += 1
							gene_num = each_name+"_gene_"+str(gene_count)+"_from_"+items[0].lower()
							scaf[items[0]] += [gene_num]
							gene_identity[gene_num] = float(items[2])
							gene_perc[gene_num] = perc
							start_end[gene_num] = [start,end]
							write_content[gene_num] = gene_num+","+ar_info[-1]+","+items[2]+\
							","+items[10]+","+ar_info[4]+","+ar_info[1]+","+ar_gene_allele+\
							","+drug_class+","+ar_mechanism+","+items[0]+","+start+","+end+"\n"

			#write the annotation content into the annotation file
			for key in write_content.keys():
				f_ar.write(write_content[key])
			print(blast_file_name+"	parse finished!!!")
			f_ar.close()

		#read the databases informations and store into a list
		annotation_info = fh.csv_reader(self.install_dir+"databases/ar_parse")
		for file in self.blast_files:
			annotation(file)

	def ar_write(self):

		def each_write(csvfilename):
			ar_genes = {}
			start_end = {}
			#get the old cds file name of genome
			old_name = self.in_dir+genome_id+".fna"
			#create the new cds file name of genome
			new_name = self.arg_dir+genome_id+"_ar"+".fna"
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
				#get the information and the location of the gene
				ar_genes[row[0]] = " ["+row[1]+"|identity:"+row[2]+"|"+row[4]+"|"+row[5]+"|position:"+row[10]+"-"+row[11]+"]"
				start_end[row[0]] = [row[9],row[10],row[11]]
			for gene in ar_genes.keys():
				for i in range(len(old_content)):
					#get the sequence of the gene according to the location in related scaffold
					if ">"+start_end[gene][0] == old_content[i].strip("\n"):
						new_content.append(">"+gene+ar_genes[gene]+"\n")
						new_content.append(old_content[i+1][int(start_end[gene][1])-1:int(start_end[gene][2])]+"\n")
			for line in new_content:
				f_new.write(line)
			f_new.close()

		#for each genome, invoke fuction each_write
		for each in self.annotation_files:
			genome_id = each.split("_ar.csv")[0]
			each_write(each)
			print("finish writing ar sequences of "+genome_id)


def main(install_directory,in_directory,out_directory):
	blastn = Blastn(install_directory,in_directory,out_directory)
	blastn.blankdel()
	blastn.makedb()
	blastn.blastn()
	pbr = ParseBlastnResult(blastn.install_dir,blastn.in_dir,blastn.out_dir)
	pbr.parse()
	pbr.ar_write()
	return blastn.out_dir

if __name__ == '__main__':
	main("","","")
