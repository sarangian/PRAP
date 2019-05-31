import os,csv,sys
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import FileHandle as fh
from sklearn.ensemble import RandomForestClassifier

class Matrix:
	#attributes of class Matrix
	#the directory of software, csv annotation files and pangenome analysis files
	install_dir = ""
	csv_file_dir = ""
	analysis_dir = ""
	#parameters for the first picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi = 200
	fontsize = 15
	glfs = 15
	fonttype = "Times New Roman"
	page_length = 20
	page_width = 15
	#the cluster parameters of the clustermap
	phenotype = "True"
	cluster_method = "average"
	row_cluster = "True"
	column_cluster = "True"
	bool_row = True
	bool_col = True
	pheno_choose = True
	column_name = "detail"

	annotation_files = []
	pheno = []
	pheno_count = {}


	def __init__(self,install_directory,csvfile_directory):
		#get the directory of software, csv annotation files
		self.install_dir = fh.cwd_get(install_directory)
		self.csv_file_dir = fh.inpd_get(csvfile_directory,file_categories="annotation csv files")
		#add all annotation files name into a list with extension
		self.annotation_files = fh.filename_get(self.csv_file_dir,"_ar.csv",showext=True)
		if len(self.annotation_files) == 0:
			print("No annotation file was found, stop analyzing...")
			sys.exit(0)
		#create the directory of pangenome analysis files
		self.analysis_dir = self.csv_file_dir+"analysis/"
		fh.dir_add(self.analysis_dir)
		#get the parameters for the clustermap
		self.dpi = int(fh.setting_reader(self.install_dir,"dpi5"))
		self.fontsize = int(fh.setting_reader(self.install_dir,"fontsize5"))
		self.fonttype = fh.setting_reader(self.install_dir,"fonttype5")
		self.glfs = fh.setting_reader(self.install_dir,"genomelabel_fsize5")
		if self.glfs == "auto":
			self.glfs = int(300/len(self.annotation_files))
			if self.glfs >= 20:
				self.glfs = 20
		else:
			self.glfs = int(self.glfs)
		self.page_length = int(fh.setting_reader(self.install_dir,"page_length5"))
		self.page_width = int(fh.setting_reader(self.install_dir,"page_width5"))
		self.cluster_method = fh.setting_reader(self.install_dir,"cluster_method5")
		self.row_cluster = fh.setting_reader(self.install_dir,"row_cluster5")
		self.column_cluster = fh.setting_reader(self.install_dir,"column_cluster5")
		self.phenotype = fh.setting_reader(self.install_dir,"show_phenotype")
		self.column_name = fh.setting_reader(self.install_dir,"column_name")
		#change the strings into bool value
		if self.row_cluster == "True":
			self.bool_row = True
		else:
			self.bool_row = False
		if self.column_cluster == "True":
			self.bool_col = True
		else:
			self.bool_col = False
		if self.phenotype == "True":
			self.pheno_choose = True
		else:
			self.pheno_choose = False

	#generate matrix file of each antibiotic in ar_phenotype.csv
	def matrix_file(self):
		#read the content in ar_phenotype.csv and return a list
		pheno_plots = fh.csv_reader(self.csv_file_dir+"res_phenotype.csv")
		#get all the phenotypes from the list
		self.pheno = [item for item in pheno_plots[0][1:]]

		#if phenotype of strain in plots
		#store the content into a dictionary
		if len(pheno_plots) > 1:
			for column in range(1,len(pheno_plots[0])):
				for row in range(1,len(pheno_plots)):
					if len(pheno_plots[row][column]) != 0:
						self.pheno_count[pheno_plots[row][0]+pheno_plots[0][column]] \
						= pheno_plots[row][column]

		#for each antibiotic in ar_phenotype.csv
		for drug_name in self.pheno:
			#the list to store antibiotic associated genes
			relate_genes = []
			write_content = {}

			print("analysis "+drug_name)
			f_matrix = open(self.analysis_dir+"3_"+drug_name+"_matrix.csv","w")

			#get all the associated genes from annotation files
			for each in self.annotation_files:
				plots = fh.csv_reader(self.csv_file_dir+each)
				for row in plots:
					if drug_name in row[7]:
						if self.column_name == "detail":
							if row[1] not in relate_genes:
								relate_genes.append(row[1])
						elif self.column_name == "allele":
							if row[6] not in relate_genes:
								relate_genes.append(row[6])
			#the title of matrix file
			write_content["title"] = "name,"
			for gene in relate_genes:
				write_content["title"] += gene+","
			write_content["title"] += "PHENOTYPE\n"
			f_matrix.write(write_content["title"])
			#get associated genes in each genome
			for each in self.annotation_files:
				#get the name of the genome
				name = each.split("_ar")[0]
				write_content[name] = name+","
				plots = fh.csv_reader(self.csv_file_dir+each)
				if self.column_name == "detail":
					column = [row[1] for row in plots]
				elif self.column_name == "allele":
					column = [row[6] for row in plots]
				#determine whether each associated gene in the genome
				for gene in relate_genes:
					if gene in column:
						row_num = column.index(gene)
						write_content[name] += str(float(plots[row_num][2])*0.01)+","
					#use "0" to represent not exist
					else:
						write_content[name] += "0,"

			#add the phenotype of relative antibiotic of the genome to the end of each line
			#if phenotype not in dictionary, use "N/A" to replace
			for key in write_content.keys():
				if key+drug_name in self.pheno_count.keys():
					write_content[key] += self.pheno_count[key+drug_name]+"\n"
				else:
					write_content[key] += "N/A\n"
			#write the content to matrix file
			for key in write_content.keys():
				if key != "title":
					f_matrix.write(write_content[key])

	#draw clustermap according to matrix files
	def matrix_picture(self):
		#change the value into float format
		def float_value(valuelist):
			for a in range(len(valuelist)):
				for b in range(len(valuelist[a])):
					valuelist[a][b] = float(valuelist[a][b])
			return valuelist

		for drug_name in self.pheno:
			try:
				print("start drawing "+drug_name+".png...")
				plots = fh.csv_reader(self.analysis_dir+"3_"+drug_name+"_matrix.csv")
				df_pheno = {}
				genome_id = []
				value = []
				#get all associated genes from matrix file
				ar_gene = plots[0][1:-1]

				#if phenotype is shown
				#only show the genome with phenotypes
				if self.pheno_choose:
					for row in plots[1:]:
						genome_id.append(row[0])
						if "N/A" not in row[-1]:
							df_pheno[row[0]] = int(row[-1])
						else:
							df_pheno[row[0]] = 0
						value.append(row[1:-1])
				else:
					for row in plots[1:]:
						genome_id.append(row[0])
						value.append(row[1:-1])
				
				value = float_value(value)
				#determine the parameters of the picture
				try:
					matplotlib.rcParams['font.sans-serif'] = self.fonttype
				except:
					print("fonttype not found!")
				matplotlib.rcParams['font.size'] = self.fontsize
				#create the picture
				if len(ar_gene) == 0:
					print("No gene/allele in "+drug_name+", drawing pirture skipped...")
					continue
				elif len(ar_gene) == 1:
					print("Only one gene/allele in "+drug_name+", drawing pirture skipped...")
					continue
				else:
					plt.subplots(figsize = (self.page_length,self.page_width))
					#create the dataframe of the clustermap
					df = pd.DataFrame(value,index=genome_id, columns=ar_gene)
					#determine the color of the map
					cmap = sns.cubehelix_palette(start = 0, dark = 0.6, light = 0.9, as_cmap = True)
					if self.pheno_choose:
						#determine the color of the phenotype
						colors = ['white', '#C3FDB8', '#B5EAAA', '#64E986', '#54C571','#4AA02C', \
						'#347C17', '#347235', '#25383C', '#254117','black']
						#determine the color of each phenotype according to the amount
						for key in df_pheno.keys():
							df_pheno[key] = colors[df_pheno[key]]
						phenotype = pd.Series(df_pheno)
						#draw the clustermap with phenotype
						f = sns.clustermap(df, cmap = cmap, row_colors = phenotype, linewidth = 0.5, \
							row_cluster = self.bool_row, col_cluster = self.bool_col, method = self.cluster_method, yticklabels = 1)
					#draw the clustermap without phenotype
					else:
						f = sns.clustermap(df, cmap = cmap, linewidth = 0.5, row_cluster = self.bool_row, \
							col_cluster = self.bool_col, method = self.cluster_method, yticklabels = 1)
					#set the format of xticklabels
					# ylabel = [i+0.5 for i in range(len(self.annotation_files))]
					# f.ax_heatmap.set_yticks(ylabel)
					plt.setp(f.ax_heatmap.get_xticklabels(), rotation=40, rotation_mode='anchor',ha="right",va="top")
					f.ax_heatmap.set_yticklabels(f.ax_heatmap.get_yticklabels(),fontsize=self.glfs)
					plt.setp(f.ax_heatmap.yaxis.get_ticklines(), markeredgewidth=10/len(self.annotation_files))
					#save the figure according to parameters
					f.savefig(self.analysis_dir+"3_"+drug_name+"_matrix.png", bbox_inches='tight', format="png", dpi=self.dpi)
					plt.close()
			except:
				print("wrong with "+drug_name+"!!!")
				pass

	def gene_contribution(self):
		if self.pheno_choose:
			for drug_name in self.pheno:
				try:
					f_cont = open(self.analysis_dir+"3_"+drug_name+"_contribution.txt","w")
					value = {}
					with open(self.csv_file_dir+"res_phenotype.csv","r") as csvfile:
						plots = list(csv.reader(csvfile, delimiter=","))
						print("start analyzing gene contribution of "+drug_name+"...")
						for row in plots[1:]:
							value[row[0]] = str(row[self.pheno.index(drug_name)+1])

					allele = []
					for name in self.annotation_files:
						with open(self.csv_file_dir+name,"r") as csvfile:
							plots = list(csv.reader(csvfile, delimiter=","))
							for row in plots[1:]:
								if row[6] not in allele:
									allele.append(row[6])
					X = []
					y = []
					for name in self.annotation_files:
						X.append([])
						y.append(value[name.split("_ar.csv")[0]])
						with open(self.csv_file_dir+name,"r") as csvfile:
							plots = list(csv.reader(csvfile, delimiter=","))
							this_allele = [row[6] for row in plots[1:]]
							for item in allele:
								X[-1].append(this_allele.count(item))
					X_train = np.array(X)
					y_train = np.array(y)

					rnd_clf = RandomForestClassifier(n_estimators=500,max_leaf_nodes=16,random_state=1)
					rnd_clf.fit(X_train,y_train)

					contribution = zip(allele,rnd_clf.feature_importances_)
					contribution = sorted(contribution,key=lambda x:x[1],reverse=True)
					print(contribution[0][0],"contribution:"+str(round(contribution[0][1]*100,2))+"%")
					for name,score in contribution:
						f_cont.write(name+"	"+str(score)+"\n")
					f_cont.close()
				except:
					print("wrong with "+drug_name+"!!!")
					pass

def main(install_directory,csvfile_directory):
	mat = Matrix(install_directory,csvfile_directory)
	if os.path.exists(mat.csv_file_dir+"res_phenotype.csv"):
		mat.matrix_file()
		if len(mat.annotation_files) == 1:
			print("Only one genome found, drawing picture skipped...")
		else:
			mat.matrix_picture()
			mat.gene_contribution()
	else:
		print("res_phenotype.csv not found, skip matrix analyze...")

if __name__ == '__main__':
	main("","")