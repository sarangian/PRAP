import FileHandle as fh
from Bio.Seq import Seq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import time
import os,sys

class Kmer:

	#attributes of class Kmer
	#including software install directory, fastq files directory
	install_dir = ""
	fastq_file_dir = ""
	kmer_dir = ""
	out_dir = ""
	fastq_files = []
	fastq_group = {}
	kmer_files = []
	annotation_files = []
	db = "res_nucl"
	k = 35
	kernel = 2
	depth = 50
	coverage = 0.95
	db_seq = {}
	db_kmer = {}
	ar_count = {}
	ar_allele = {}
	drug_class = {}

	def __init__(self,install_directory,fastq_directory,out_directory):

		def SampleGroup(files):
			for file in files:
				group = ""
				if file.endswith(".1.fastq"):
					group = file.split(".1.fastq")[0]
				elif file.endswith(".2.fastq"):
					group = file.split(".2.fastq")[0]
				else:
					if file.endswith(".fastq"):
						group = file.split(".fastq")[0]

				if group not in self.fastq_group:
					self.fastq_group[group] = [file]
				else:
					self.fastq_group[group].append(file)

		#get the install directory of the software
		self.install_dir = fh.cwd_get(install_directory)
		#get the fasta files directory
		self.fastq_file_dir = fh.inpd_get(fastq_directory,file_categories="fastq files")
		self.out_dir = fh.inpd_get(out_directory,file_categories="ouput files")
		#get fastq files
		self.fastq_files = fh.filename_get(self.fastq_file_dir,".fastq",showext = True)
		if len(self.fastq_files) == 0:
			print("No fastq file was found, stop analyzing...")
			sys.exit(0)
		self.kmer_dir = self.out_dir+"kmer/"
		fh.dir_add(self.kmer_dir)
		self.annotation_files = fh.filename_get(self.out_dir,"_ar.csv",showext = True)
		SampleGroup(self.fastq_files)
		#get kmer
		self.k = int(fh.setting_reader(self.install_dir,"k_value"))
		if not bool(self.k%2):
			print("k value must be an odd number!")
			sys.exit(0)
		self.kernel = int(fh.setting_reader(self.install_dir,"kernel"))
		self.depth = int(fh.setting_reader(self.install_dir,"depth"))
		self.area = int(fh.setting_reader(self.install_dir,"area_score"))
		self.coverage = float(fh.setting_reader(self.install_dir,"kmer_coverage"))

		localtime = time.asctime(time.localtime(time.time()))
		print("kmer start time: ",localtime)

	def KmerDB(self):

		def AlleleGet():
			plots = fh.csv_reader(self.install_dir+"databases/res_parse")
			for row in plots:
				self.ar_allele[row[0]] = row[2]
				self.drug_class[row[0]] = row[1]

		def DBget(db_name):
			print("reading database: %s \n" % db_name)
			with open(db_name,"r") as f:
				content = f.readlines()
				content.append("")
				indices = [i for i in range(len(content)) if ">" in content[i]]
				indices.append(len(content))
				for i in range(len(indices)-1):
					seq_name = content[indices[i]].split(">")[1].split("\n")[0]
					self.db_seq[seq_name] = ""
					for l in range(indices[i]+1,indices[i+1]):
						self.db_seq[seq_name] += content[l].strip("\n")

		def DBform():
			print("forming kmer (k=%d) database, please wait..." % self.k)
			for key in self.db_seq.keys():
				seq = Seq(self.db_seq[key])
				reco_seq = seq.reverse_complement()
				kmer_count = len(self.db_seq[key])-self.k+1
				for i in range(kmer_count):
					if str(seq[i:i+self.k]) not in self.db_kmer:
						self.db_kmer[str(seq[i:i+self.k])] = [[key,i]]
					else:
						self.db_kmer[str(seq[i:i+self.k])].append([key,i])
					if str(reco_seq[i:i+self.k]) not in self.db_kmer:
						self.db_kmer[str(reco_seq[i:i+self.k])] = [[key,i]]
					else:
						self.db_kmer[str(reco_seq[i:i+self.k])].append([key,i])
			print("totally %d database reads (k=%d) are formed...\n" % (len(self.db_kmer),self.k))

		AlleleGet()
		DBget(self.install_dir+"databases/"+self.db)
		DBform()

	def KmerReads(self):

		def ReadsGet(fastq_file):

			def UseReads(seq):

				length = len(seq)
				radius = self.k//2

				if self.kernel == 1:
					if len(seq) >= self.k:
						mid = length//2
						mid_seq = ""
						if bool(self.k%2):
							mid_seq = seq[mid-radius:mid+radius+1]

						if mid_seq in self.db_kmer:
							return True
						else:
							return False
					else:
						print("k value should be smaller than the reads' length!!!")
						sys.exit(0)

				elif self.kernel == 2:
					if float(len(seq)) >= 1.5*self.k:
						mid_1 = length//3
						mid_2 = length*2//3
						mid_seq_1 = ""
						mid_seq_2 = ""

						if bool(self.k%2):
							mid_seq_1 = seq[mid_1-radius:mid_1+radius+1]
							mid_seq_2 = seq[mid_2-radius:mid_2+radius+1]

						if mid_seq_1 in self.db_kmer or mid_seq_2 in self.db_kmer:
							return True
						else:
							return False
					else:
						print("k value should be smaller than two-thirds of the reads' length")
						sys.exit(0)

				elif self.kernel == 3:
					if len(seq) >= 2*self.k:
						mid = length//2
						mid_1 = length//4
						mid_2 = length*3//4
						mid_seq = ""
						mid_seq_1 = ""
						mid_seq_2 = ""

						if bool(self.k%2):
							mid_seq = seq[mid-radius:mid+radius+1]
							mid_seq_1 = seq[mid_1-radius:mid_1+radius+1]
							mid_seq_2 = seq[mid_2-radius:mid_2+radius+1]

						if mid_seq_1 in self.db_kmer or mid_seq_2 in self.db_kmer or mid_seq in self.db_kmer:
							return True
						else:
							return False
					else:
						print("k value should be smaller than a half of the reads' length")
						sys.exit(0)

				else:
					print("the value of kernel should be 1, 2 or 3, please select the write number")
					sys.exit(0)

			print("kmerizing reads in %s (k=%d)..." % (fastq_file, self.k))
			with open(self.fastq_file_dir+fastq_file,"r") as f:
				marker = 0
				while True:
					line = f.readline()
					if not line:
						break
					if marker == 1:
						sequence = line.strip("\n")
						use_read = UseReads(sequence)

						if use_read:
							for i in range(len(sequence)-self.k+1):
								if sequence[i:i+self.k] in self.db_kmer:
									for gene in self.db_kmer[sequence[i:i+self.k]]:
										self.ar_count[self.ar_allele\
										[gene[0].split("_")[0]]][gene[0]][gene[1]] += 1
						marker = 0
					if "@" in line:
						marker = 1
		try:
			plt.rcParams['font.family'] = "Times New Roman"
		except:
			print("fonttype not found!")
		for group in self.fastq_group:

			for key in self.db_seq.keys():
				gene = key.split("_")[0]
				if self.ar_allele[gene] not in self.ar_count:
					self.ar_count[self.ar_allele[gene]] = {}
					self.ar_count[self.ar_allele[gene]][key] = {}
				else:
					self.ar_count[self.ar_allele[gene]][key] = {}
				for i in range(len(self.db_seq[key])-self.k+1):
					self.ar_count[self.ar_allele[gene]][key][i] = 0
					
			print("analyzing group %s..." % group)
			predicted_gene = {}
			for each in self.fastq_group[group]:
				ReadsGet(each)

			print("writing results...")
			f_ar = open(self.kmer_dir+group+"_"+str(self.k)+"mer_countVS"+self.db+".csv","w")
			self.kmer_files.append(group+"_"+str(self.k)+"mer_countVS"+self.db+".csv")
			for allele in self.ar_count:
				gene_total = {}
				gene_perc = {}
				best_gene = ""
				y_max = 0
				for gene in self.ar_count[allele]:
					total = 0
					freq = 0
					former_loc = ""
					fall = False
					for loc in self.ar_count[allele][gene]:
						try:
							minus = self.ar_count[allele][gene][loc]-self.ar_count[allele][gene][former_loc]
							if abs(minus) >= 50:
								fall = True
						except:
							pass
						former_loc = loc
						if self.ar_count[allele][gene][loc] >= self.depth:
							if self.ar_count[allele][gene][loc] > y_max:
								y_max = self.ar_count[allele][gene][loc]
							total += self.ar_count[allele][gene][loc]
							freq += 1
					total = total/len(self.ar_count[allele][gene])
					if total >= self.area:
						perc = (freq/len(self.ar_count[allele][gene]))*100
						if perc >= self.coverage*100:
							if not fall:
								best_gene = gene
							gene_total[gene] = total
							gene_perc[gene] = perc
							f_ar.write(gene+","+str(total))
							f_ar.write(","+str(perc)+"%\n")
				if len(gene_total) >= 1:
					f,ax = plt.subplots(1,1)
					max_gene = sorted(gene_perc.items(),key = lambda x:x[1],reverse=True)
					max_perc = max_gene[0][1]
					most_prob_gene = max_gene[0][0]
					for each in max_gene:
						if gene_perc[each[0]] == max_perc and gene_total[each[0]] > gene_total[most_prob_gene]:
							most_prob_gene = each[0]
					most_index = max_gene.index((most_prob_gene,max_perc))
					if most_index != 0:
						max_gene[0],max_gene[most_index] = max_gene[most_index],max_gene[0]
					genes = [gene[0] for gene in max_gene]
					if best_gene != "":
						best_index = genes.index(best_gene)
						if best_index != 0:
							genes[0],genes[best_index] = genes[best_index],genes[0]
							most_prob_gene = best_gene
					predicted_gene[most_prob_gene] = \
					[gene_total[most_prob_gene],gene_perc[most_prob_gene],allele]

					if len(genes) > 4:
						genes = genes[0:4]
					for gene in genes:
						x_list = []
						y_list = []
						for loc in self.ar_count[allele][gene]:
							x_list.append(loc)
							y_list.append(self.ar_count[allele][gene][loc])
						if genes.index(gene) == 0:
							plt.plot(x_list,y_list,label=gene,linewidth=1.0,color="black")
						else:
							plt.plot(x_list,y_list,label=gene,linewidth=0.5)

						plt.legend(loc='best')
						plt.ylim((0,y_max+50))
						ax.set_title(str(allele))
						ax.set_xlabel("gene length (bp)")
						ax.set_ylabel(str(self.k)+"mer reads count")
						if not os.path.exists(self.kmer_dir+group+"/"):
							os.mkdir(self.kmer_dir+group+"/")
						if not os.path.exists(self.kmer_dir+group+"/"+self.db+"_"+str(self.k)+"/"):
							os.mkdir(self.kmer_dir+group+"/"+self.db+"_"+str(self.k)+"/")
						f.savefig(self.kmer_dir+group+"/"+self.db+"_"+str(self.k)+"/"+allele+"_"+str(self.k)+"mer.png", \
							bbox_inches='tight', format="png", dpi=200)
					plt.close()
			f_ar.close()

			f_pred = open(self.out_dir+group+"_ar.csv","w")
			f_pred.write("#gene_num"+","+"gene_name"+","+"coverage"+","+"area_score"\
				+","+"ARO_num"+","+"accession_num"+","+"ar_gene_allele"+","+"drug_class"\
				+","+"ar_mechanism"+"\n")
			gene_count = 0
			for gene in predicted_gene:
				gene_num = group+"_AMRgene_"+str(gene_count)
				gene_count += 1
				accession_num = gene.split("_")[-1] if "NC" not in gene else "NC_"+gene.split("_")[-1]
				f_pred.write(gene_num+","+gene+","+str(predicted_gene[gene][1])+","+str(predicted_gene[gene][0])+\
					",/,"+accession_num+","+predicted_gene[gene][2]+","+self.drug_class[gene.split("_")[0]]+",/\n")
			f_pred.close()
			print("finish analyzing group %s...\n" % group)

def main(install_directory,fastq_directory,out_directory):
	kmer = Kmer(install_directory,fastq_directory,out_directory)
	if len(kmer.annotation_files) == len(kmer.fastq_group):
		choose = input("kmer files exist, do you want to update them?(Y/N):").upper()
		if choose == "Y":
			kmer.KmerDB()
			kmer.KmerReads()
	else:
		kmer.KmerDB()
		kmer.KmerReads()
	localtime = time.asctime(time.localtime(time.time()))
	print("kmer end time: ",localtime)
	return kmer.out_dir

if __name__ == '__main__':
	main("","","")