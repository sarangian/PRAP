import os,csv,sys
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import FileHandle as fh
	
#this part is used to summarize and analyze the pan resistome
class Pansum:

	#attributes of class Pansum
	#the directory of software, csv annotation files and analysis files
	install_dir = ""
	csv_file_dir = ""
	analysis_dir = ""
	drug_class  = []
	mech_class = []
	drug_count = {}
	mech_count = {}
	#parameters for summary picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi = 200
	fontsize = 20
	glfs = 20
	fonttype = "Times New Roman"
	page_length = 20
	page_width = 15
	#parameters for the first picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi1 = 200
	fontsize1 = 20
	glfs1 = 20
	fonttype1 = "Times New Roman"
	page_length1 = 20
	page_width1 = 15
	#the cluster parameters of the clustermap
	cluster_method = "average"
	row_cluster = "True"
	column_cluster = "True"
	bool_row = True
	bool_col = True
	#parameters for the second picture
	dpi2 = 200
	fontsize2 = 15
	glfs2 = 15
	fonttype2 = "Times New Roman"
	page_length2 = 20
	page_width2 = 15
	dotsize = 30

	def __init__(self,install_directory,csvfile_directory):
		#get the directory of software, csv annotation files
		self.install_dir = fh.cwd_get(install_directory)
		self.csv_file_dir = fh.inpd_get(csvfile_directory,file_categories="annotation csv files")
		#add all annotation files name into a list with extension
		self.annotation_files = fh.filename_get(self.csv_file_dir,"_ar.csv",showext=True)
		if len(self.annotation_files) == 0:
			print("No annotation file was found, stop analyzing...")
			sys.exit(0)
		#create the directory of analysis files
		self.analysis_dir = self.csv_file_dir+"analysis/"
		fh.dir_add(self.analysis_dir)
		#get the parameters for the picture from settings.txt (1)ar distribution picture
		self.dpi = int(fh.setting_reader(self.install_dir,"dpi0"))
		self.fontsize = int(fh.setting_reader(self.install_dir,"fontsize0"))
		self.glfs = fh.setting_reader(self.install_dir,"genomelabel_fsize0")
		if self.glfs == "auto":
			self.glfs = int(400/len(self.annotation_files))
			if self.glfs >= 30:
				self.glfs = 30
		else:
			self.glfs = int(self.glfs)
		self.fonttype = fh.setting_reader(self.install_dir,"fonttype0")
		self.page_length = int(fh.setting_reader(self.install_dir,"page_length0"))
		self.page_width = int(fh.setting_reader(self.install_dir,"page_width0"))

		#get the parameters for the clustermap
		self.dpi1 = int(fh.setting_reader(self.install_dir,"dpi3"))
		self.fontsize1 = int(fh.setting_reader(self.install_dir,"fontsize3"))
		self.glfs1 = fh.setting_reader(self.install_dir,"genomelabel_fsize3")
		if self.glfs1 == "auto":
			self.glfs1 = int(400/len(self.annotation_files))
			if self.glfs1 >= 30:
				self.glfs1 = 30
		else:
			self.glfs1 = int(self.glfs1)
		self.fonttype1 = fh.setting_reader(self.install_dir,"fonttype3")
		self.page_length1 = int(fh.setting_reader(self.install_dir,"page_length3"))
		self.page_width1 = int(fh.setting_reader(self.install_dir,"page_width3"))
		self.cluster_method = fh.setting_reader(self.install_dir,"cluster_method3")
		self.row_cluster = fh.setting_reader(self.install_dir,"row_cluster3")
		self.column_cluster = fh.setting_reader(self.install_dir,"column_cluster3")
		#change the strings into bool value
		if self.row_cluster == "True":
			self.bool_row = True
		else:
			self.bool_row = False
		if self.column_cluster == "True":
			self.bool_col = True
		else:
			self.bool_col = False

		#get the parameters of the correlation map
		self.dpi2 = int(fh.setting_reader(self.install_dir,"dpi4"))
		self.dotsize = int(fh.setting_reader(self.install_dir,"dot_size4"))
		self.fontsize2 = int(fh.setting_reader(self.install_dir,"fontsize4"))
		self.glfs2 = fh.setting_reader(self.install_dir,"genomelabel_fsize4")
		if self.glfs2 == "auto":
			self.glfs2 = int(400/len(self.annotation_files))
			if self.glfs2 >= 30:
				self.glfs2 = 30
		else:
			self.glfs2 = int(self.glfs2)
		self.fonttype2 = fh.setting_reader(self.install_dir,"fonttype4")
		self.page_length2 = int(fh.setting_reader(self.install_dir,"page_length4"))
		self.page_width2 = int(fh.setting_reader(self.install_dir,"page_width4"))

	#generate a summary file of pan antibiotoc resistance genome
	def pansum_file(self):
		print("start summary...")
		#create the file of summary of drug classes and mechanism classes
		f_sum_ar = open(self.analysis_dir+"1_class_summary.csv","w")
		f_sum_mech = open(self.analysis_dir+"1_mech_summary.csv","w")
		#the dictionary to store contents which will be written
		write_content_ar = {}
		write_content_mech = {}

		#a function to get all the drug class from the annotation files
		#return a list of drug class
		def title_get(title_list,plot,p):
			for row in plot[1:]:
				#if more than one drug in an item
				if ";" in row[p]:
					items = row[p].split(";")
					for item in items:
						if item not in title_list:
							title_list.append(item)
				else:
					if row[p] not in title_list:
						title_list.append(row[p])
			title_list.sort()
			return title_list

		#a function to count all the items given in the parameters:title
		#return a dictionaty which keys are items in title 
		#and values related to amount of the item
		def count_get(title,plot,p):
			count_dict = {}
			#give each item a initial value of 0
			for key in title:
				count_dict[key] = 0

			for row in plot[1:]:
				#if more than one drug in an item
				if ";" in row[p]:
					items = row[p].split(";")
					for item in items:
						count_dict[item] += 1
				else:
					count_dict[row[p]] += 1
			return count_dict

		#get all the drug classes and mechanism classes for annotation fils
		for each in self.annotation_files:
			plots = fh.csv_reader(self.csv_file_dir+each)
			self.drug_class = title_get(self.drug_class,plots,7)
			self.mech_class = title_get(self.mech_class,plots,8)
		#write the title of the summary file
		write_content_ar["title"] = "name,"
		write_content_mech["title"] = "name,"
		for each in self.drug_class:
			write_content_ar["title"] += each+","
		for each in self.mech_class:
			write_content_mech["title"] += each+","
		write_content_ar["title"] = write_content_ar["title"].strip(",")
		write_content_mech["title"] = write_content_mech["title"].strip(",")
		write_content_ar["title"] += "\n"
		write_content_mech["title"] += "\n"
		f_sum_ar.write(write_content_ar["title"])
		f_sum_mech.write(write_content_mech["title"])

		#count the number of drug classes and mechanism classes
		for each in self.annotation_files:
			print("counting "+each+"...")
			name = each.split("_ar")[0]
			plots = fh.csv_reader(self.csv_file_dir+each)
			self.drug_count = count_get(self.drug_class,plots,7)
			self.mech_count = count_get(self.mech_class,plots,8)
			#write the amount according to the genome id and drug class
			write_content_ar[name] = name+","
			write_content_mech[name] = name+","
			for each_drug in self.drug_class:
				write_content_ar[name] += str(self.drug_count[each_drug])+","
			for each_mech in self.mech_class:
				write_content_mech[name] += str(self.mech_count[each_mech])+","
			write_content_ar[name] = write_content_ar[name].strip(",")
			write_content_mech[name] = write_content_mech[name].strip(",")
			write_content_ar[name] += "\n"
			write_content_mech[name] += "\n"
		#write all the contents into the summary file
		for key in write_content_ar.keys():
			if key != "title":
				f_sum_ar.write(write_content_ar[key])
		for key in write_content_mech.keys():
			if key != "title":
				f_sum_mech.write(write_content_mech[key])
		f_sum_ar.close()
		f_sum_mech.close()

	#to draw the bar graph of the pan ar genome
	def pansum_picture(self):
		print("drawing bar graph for pan-resistome...")
		value = {}
		#get the datasets from the summary file
		class_plots = fh.csv_reader(self.analysis_dir+"1_class_summary.csv")
		#change the string format to integer
		for i in range(1,len(class_plots[0])):
			value[class_plots[0][i]] = [int(row[i]) for row in class_plots[1:]]
		#get the x_axis of the bar graph
		x_axis = [row[0] for row in class_plots[1:]]
		#the list to store different classes
		y_hue = []
		#determine the fontsize and fonttype according to settings
		try:
			plt.rcParams['font.family'] = self.fonttype
		except:
			print("fonttype not found!")
		plt.rcParams['font.size'] = self.fontsize
		#determine the figure size according to settings
		f,ax = plt.subplots(figsize = (self.page_length,self.page_width))
		#determine the xticklabels and rotate the labels
		ax.set_xticklabels(x_axis, rotation=40, fontsize=self.glfs,\
		 rotation_mode='anchor', ha="right", va="top")
		#create the bottom of the bar
		bottom_class = [0 for i in range(len(x_axis))]
		#the list of the bar already drawn
		axes = []
		for key in value:
			if key != "name":
				#draw the bar of current drug
				p = ax.bar(x_axis,value[key],bottom=bottom_class)
				y_hue.append(key)
				axes.append(p)
				#add the bottom of the bar
				for i in range(len(bottom_class)):
					bottom_class[i] += value[key][i]
		#create the legend of the graph
		plt.legend(axes[::-1],y_hue[::-1],loc=3,bbox_to_anchor = (1.02, 0))
		#set name of x axis
		ax.set_xlabel('Genome Id')
		#set name of y axis
		ax.set_ylabel('Antibiotic Resistance Genes Count')
		#save the figure
		ax.set_title('Overview of Pan-resistome')
		f.savefig(self.analysis_dir+"1_class_summary.png",\
		 bbox_inches='tight', format="png", dpi=self.dpi)
		plt.close()

	#draw pictures of the ar genes
	def pan_picture(self):
		print("drawing clustermap for pan-resistome...")
		#get the dataset from the summary file
		value = []
		plots = fh.csv_reader(self.analysis_dir+"1_class_summary.csv")
		#get the genome ids from the summary file
		genome_id = [plots[row][0] for row in range(1,len(plots))]
		#get the dataset of the value from the plots
		#each list in the lists of value represent each row of plots
		#each item in each list represent elements in each row
		for row in range(1,len(plots)):
			value.append([])
			for column in range(1,len(plots[0])):
				value[row-1].append(int(plots[row][column]))
		#transfer the list into array format using np.array()
		genome_id = np.array(genome_id)
		self.drug_class = np.array(self.drug_class)
		value = np.array(value)
		subplots_num = len(genome_id)

		#draw clustermap
		#set the parameters of the picture
		try:
			plt.rcParams['font.family'] = self.fonttype1
		except:
			print("fonttype not found!")
		plt.rcParams['font.size'] = self.fontsize1
		#the dataframe of the clustermap
		df = pd.DataFrame(value,index=genome_id, columns=self.drug_class)
		plt.subplots(figsize = (self.page_length1,self.page_width1))
		#the color of the clustermap
		cmap = sns.cubehelix_palette(start = 0, rot = -0.1,\
		 dark = 0.3, light = 1.0, as_cmap = True)
		#draw the clustermap
		f = sns.clustermap(df, cmap = cmap, linewidth = 0.5,\
		 row_cluster = self.bool_row, yticklabels = 1, \
		 col_cluster = self.bool_col, method = self.cluster_method)
		#set the format of the xticklables and rotate them
		plt.setp(f.ax_heatmap.get_xticklabels(), rotation=40,\
		 rotation_mode='anchor',ha="right",va="top")
		plt.setp(f.ax_heatmap.get_yticklabels(), \
			fontsize=self.glfs1)
		plt.setp(f.ax_heatmap.yaxis.get_ticklines(), \
			markeredgewidth=10/len(self.annotation_files))
		#save the clustermap
		f.savefig(self.analysis_dir+"1_ar_cluster.png",\
		 bbox_inches='tight', format="png", dpi=self.dpi1)
		plt.close()
		print("finish drawing ARGs cluster...")

		if len(self.annotation_files) <= 60:
			#draw correlation map
			#set the parameters of the picture
			try:
				plt.rcParams['font.family'] = self.fonttype2
			except:
				print("fonttype not found!")
			plt.rcParams['font.size'] = self.glfs2
			#draw the a blank canvas with number of n^2 subplots
			#(n refers to the number of genomes)
			f,axs = plt.subplots(nrows = subplots_num, ncols = subplots_num,\
				sharex= True, sharey = True, figsize = (self.page_length2,self.page_width2))
			f.subplots_adjust(wspace=0, hspace=0)
			axs[0,0].invert_xaxis()

			#draw the picture in each subplot
			for i in range(subplots_num**2):
				#show the progress of the process
				print("%d/%d drawing No.%d subgraph of total %d subgraphs"\
				 % (i+1,subplots_num**2,i+1,subplots_num**2))
				color_index = []
				#get the position of the subplot
				row_num = i//subplots_num
				col_num = i%subplots_num
				ax = axs[row_num,col_num]
				#give different color of different dot
				for j in range(len(value[row_num])):
					if value[row_num][j] > value[col_num][j]:
						color_index.append("x>y")
					else:
						if value[row_num][j] < value[col_num][j]:
							color_index.append("x<y")
						else:
							color_index.append("x=y")
				#draw the scatterplot map
				sns.scatterplot(x = value[col_num], y = value[row_num], \
					hue = color_index, style = color_index, \
					s= int(self.dotsize), \
					ax = ax, \
					hue_order = ["x=y","x>y","x<y"], style_order = ["x>y","x=y","x<y"], \
					legend = False)
				#rotate the y labels
				if col_num == 0:
					axs[row_num,col_num].set_ylabel(genome_id[row_num], \
						rotation = 0, rotation_mode='anchor', ha="right", va="top")
				#rotate the x labels
				if row_num == len(genome_id)-1:
					axs[row_num,col_num].set_xlabel(genome_id[col_num], \
						rotation=40, rotation_mode='anchor', ha="right", va="top")
			#save the correlation map
			print("saving figure...")
			f.savefig(self.analysis_dir+"1_ar_corr.png", \
				bbox_inches='tight', format="png", dpi=self.dpi2)
			plt.close()
			print("finish drawing ARGs correlation graph...")
		elif len(self.annotation_files) > 60:
			print("more than 60 genomes were input, correlation graph analysis skipped...")

#this part is used to summarize and analyze the accessory resistome
class AccessorySummary:
	#attributes of class AccessorySummary
	#the directory of software, csv annotation files and pangenome analysis files
	install_dir = ""
	csv_file_dir = ""
	analysis_dir = ""
	#parameters for summary picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi = 200
	fontsize = 20
	glfs = 20
	fonttype = "Times New Roman"
	page_length = 20
	page_width = 15
	#parameters for the first picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi1 = 200
	fontsize1 = 20
	glfs1 = 20
	fonttype1 = "Times New Roman"
	page_length1 = 20
	page_width1 = 15
	#the cluster parameters of the clustermap
	cluster_method = "average"
	row_cluster = "True"
	column_cluster = "True"
	bool_row = True
	bool_col = True
	#parameters for the second picture
	dpi2 = 200
	fontsize2 = 15
	glfs2 = 15
	fonttype2 = "Times New Roman"
	page_length2 = 20
	page_width2 = 15
	dotsize = 30

	annotation_files = []
	ar_genes = {}
	accessory_ar_genes = {}
	drug_title = []
	drug_class = {}

	def __init__(self,install_directory,csvfile_directory):
		#get the directory of software, csv annotation files
		self.install_dir = fh.cwd_get(install_directory)
		self.csv_file_dir = fh.inpd_get(csvfile_directory,file_categories="annotation csv files")
		#create the directory of pangenome analysis files
		self.analysis_dir = self.csv_file_dir+"analysis/"
		fh.dir_add(self.analysis_dir)
		#add all annotation files name into a list with extension
		self.annotation_files = fh.filename_get(self.csv_file_dir,"_ar.csv",showext=True)

		self.dpi = int(fh.setting_reader(self.install_dir,"dpi0"))
		self.fontsize = int(fh.setting_reader(self.install_dir,"fontsize0"))
		self.glfs = fh.setting_reader(self.install_dir,"genomelabel_fsize0")
		if self.glfs == "auto":
			self.glfs = int(400/len(self.annotation_files))
			if self.glfs >= 30:
				self.glfs = 30
		else:
			self.glfs = int(self.glfs)
		self.fonttype = fh.setting_reader(self.install_dir,"fonttype0")
		self.page_length = int(fh.setting_reader(self.install_dir,"page_length0"))
		self.page_width = int(fh.setting_reader(self.install_dir,"page_width0"))

		#get the parameters for the clustermap
		self.dpi1 = int(fh.setting_reader(self.install_dir,"dpi3"))
		self.fontsize1 = int(fh.setting_reader(self.install_dir,"fontsize3"))
		self.glfs1 = fh.setting_reader(self.install_dir,"genomelabel_fsize3")
		if self.glfs1 == "auto":
			self.glfs1 = int(400/len(self.annotation_files))
			if self.glfs1 >= 30:
				self.glfs1 = 30
		else:
			self.glfs1 = int(self.glfs1)
		self.fonttype1 = fh.setting_reader(self.install_dir,"fonttype3")
		self.page_length1 = int(fh.setting_reader(self.install_dir,"page_length3"))
		self.page_width1 = int(fh.setting_reader(self.install_dir,"page_width3"))
		self.cluster_method = fh.setting_reader(self.install_dir,"cluster_method3")
		self.row_cluster = fh.setting_reader(self.install_dir,"row_cluster3")
		self.column_cluster = fh.setting_reader(self.install_dir,"column_cluster3")
		#change the strings into bool value
		if self.row_cluster == "True":
			self.bool_row = True
		else:
			self.bool_row = False
		if self.column_cluster == "True":
			self.bool_col = True
		else:
			self.bool_col = False

		#get the parameters of the correlation map
		self.dpi2 = int(fh.setting_reader(self.install_dir,"dpi4"))
		self.dotsize = int(fh.setting_reader(self.install_dir,"dot_size4"))
		self.fontsize2 = int(fh.setting_reader(self.install_dir,"fontsize4"))
		self.glfs2 = fh.setting_reader(self.install_dir,"genomelabel_fsize4")
		if self.glfs2 == "auto":
			self.glfs2 = int(400/len(self.annotation_files))
			if self.glfs2 >= 30:
				self.glfs2 = 30
		else:
			self.glfs2 = int(self.glfs2)
		self.fonttype2 = fh.setting_reader(self.install_dir,"fonttype4")
		self.page_length2 = int(fh.setting_reader(self.install_dir,"page_length4"))
		self.page_width2 = int(fh.setting_reader(self.install_dir,"page_width4"))

	#generate the accessory ar genes files
	def access_file(self):
		#get all the ar genes and their amounts
		ar_genes = {}
		for each in self.annotation_files:
			plots = fh.csv_reader(self.csv_file_dir+each)
			for row in plots[1:]:
				if row[1] not in ar_genes:
					ar_genes[row[1]] = 1
				else:
					ar_genes[row[1]] += 1

		#get all the accessory ar genes according to their amounts 
		genome_num = len(self.annotation_files)
		accessory_ar_genes = {}
		for gene in ar_genes.keys():
			#if the amount of the gene is smaller than the amount of genomes
			if ar_genes[gene] < genome_num:
				accessory_ar_genes[gene] = ar_genes[gene]

		#write the file of each accessory ar genes
		for each in self.annotation_files:
			print("analyzing accessory genes of "+each)
			name = each.split("_ar")[0]
			#create the accessory ar genes file
			f_accessory = open(self.csv_file_dir+name+"_ar_accessory.csv","w")
			#read the original ar gene files
			plots = fh.csv_reader(self.csv_file_dir+each)
			for row in plots:
				write_content = ""
				if row[1] in accessory_ar_genes.keys():
					for item in row:
						write_content += item+","
					write_content = write_content.strip(",")+"\n"
					f_accessory.write(write_content)
			f_accessory.close()

	#summarize the accessory ar genes files
	def accsum_file(self):
		print("\naccessory summary start...")
		#create the summary file of accessory ar genes
		f_accessory_summary = open(self.analysis_dir+"2_accessory_class_summary.csv","w")
		#get all the drug class
		for each in self.annotation_files:
			name = each.split("_ar")[0]
			plots = fh.csv_reader(self.csv_file_dir+name+"_ar_accessory.csv")
			for row in plots[1:]:
				#if more than one drug in the drug class
				if ";" in row[7]:
					#split the drugs into list of each drug
					drugs = row[7].split(";")
					for drug in drugs:
						if drug not in self.drug_title:
							self.drug_title.append(drug)
				else:
					if row[7] not in self.drug_title:
						self.drug_title.append(row[7])

		#write the title of the summary file
		self.drug_title.sort()
		write_line = "name,"
		for each in self.drug_title:
			write_line += each+","
		write_line = write_line.strip(",")+"\n"
		f_accessory_summary.write(write_line)

		#write the content of the summary file
		for each in self.annotation_files:
			name = each.split("_ar")[0]
			#give the initial value of 0
			for each_one in self.drug_title:
				self.drug_class[each_one] = 0
			plots = fh.csv_reader(self.csv_file_dir+name+"_ar_accessory.csv")
			#count for the amount
			for row in plots[1:]:
				#if more than one drug in the drug class
				if ";" in row[7]:
					#split the drugs into list of each drug
					drugs = row[7].split(";")
					for drug in drugs:
						self.drug_class[drug] += 1
				else:
					self.drug_class[row[7]] += 1
			#write the accessory ar genes information for each genome
			write_line = name+","
			for each in self.drug_title:
				write_line += str(self.drug_class[each])+","
			write_line = write_line.strip(",")+"\n"
			f_accessory_summary.write(write_line)
		f_accessory_summary.close()
		print("accessory summary finished...\n")

	#to draw the bar graph of the pan ar genome
	def accsum_picture(self):
		print("drawing bar graph for accessory resistome...")
		value = {}
		#get the datasets from the summary file
		class_plots = fh.csv_reader(self.analysis_dir+"2_accessory_class_summary.csv")
		#change the string format to integer
		for i in range(1,len(class_plots[0])):
			value[class_plots[0][i]] = [int(row[i]) for row in class_plots[1:]]
		#get the x_axis of the bar graph
		x_axis = [row[0] for row in class_plots[1:]]
		#the list to store different classes
		y_hue = []
		#determine the fontsize and fonttype according to settings
		try:
			plt.rcParams['font.family'] = self.fonttype
		except:
			print("fonttype not found!")
		plt.rcParams['font.size'] = self.fontsize
		#determine the figure size according to settings
		f,ax = plt.subplots(figsize = (self.page_length,self.page_width))
		#determine the xticklabels and rotate the labels
		ax.set_xticklabels(x_axis, rotation=40, fontsize=self.glfs,\
		 rotation_mode='anchor', ha="right", va="top")
		#create the bottom of the bar
		bottom_class = [0 for i in range(len(x_axis))]
		#the list of the bar already drawn
		axes = []
		for key in value:
			if key != "name":
				#draw the bar of current drug
				p = ax.bar(x_axis,value[key],bottom=bottom_class)
				y_hue.append(key)
				axes.append(p)
				#add the bottom of the bar
				for i in range(len(bottom_class)):
					bottom_class[i] += value[key][i]
		#create the legend of the graph
		plt.legend(axes[::-1],y_hue[::-1],loc=3,bbox_to_anchor = (1.02, 0))
		#set name of x axis
		ax.set_xlabel('Genome Id')
		#set name of y axis
		ax.set_ylabel('Antibiotic Resistance Genes Count')
		#save the figure
		ax.set_title('Overview of Accessory Resistomes')
		f.savefig(self.analysis_dir+"2_accessory_class_summary.png",\
		 bbox_inches='tight', format="png", dpi=self.dpi)
		plt.close()

	#draw pictures of the accessory ar genes
	def acc_picture(self):
		print("drawing clustermap for accessory resistome...")
		#get the dataset from the summary file
		value = []
		plots = fh.csv_reader(self.analysis_dir+"2_accessory_class_summary.csv")
		#get the genome ids from the summary file
		genome_id = [plots[row][0] for row in range(1,len(plots))]
		#get the dataset of the value from the plots
		#each list in the lists of value represent each row of plots
		#each item in each list represent elements in each row
		for row in range(1,len(plots)):
			value.append([])
			for column in range(1,len(plots[0])):
				value[row-1].append(int(plots[row][column]))
		#transfer the list into array format using np.array()
		genome_id = np.array(genome_id)
		self.drug_title = np.array(self.drug_title)
		value = np.array(value)
		subplots_num = len(genome_id)

		#draw clustermap
		#set the parameters of the picture
		try:
			plt.rcParams['font.family'] = self.fonttype1
		except:
			print("fonttype not found!")
		plt.rcParams['font.size'] = self.fontsize1
		#the dataframe of the clustermap
		df = pd.DataFrame(value,index=genome_id, columns=self.drug_title)
		plt.subplots(figsize = (self.page_length1,self.page_width1))
		#the color of the clustermap
		cmap = sns.cubehelix_palette(start = 0, rot = -0.1,\
		 dark = 0.3, light = 1.0, as_cmap = True)
		#draw the clustermap
		f = sns.clustermap(df, cmap = cmap, linewidth = 0.5,\
		 row_cluster = self.bool_row, yticklabels = 1, \
		 col_cluster = self.bool_col, method = self.cluster_method)
		#set the format of the xticklables and rotate them
		plt.setp(f.ax_heatmap.get_xticklabels(), rotation=40,\
		 rotation_mode='anchor',ha="right",va="top")
		plt.setp(f.ax_heatmap.get_yticklabels(), \
			fontsize=self.glfs1)
		plt.setp(f.ax_heatmap.yaxis.get_ticklines(), \
			markeredgewidth=10/len(self.annotation_files))
		#save the clustermap
		f.savefig(self.analysis_dir+"2_accessory_ar_cluster.png",\
		 bbox_inches='tight', format="png", dpi=self.dpi1)
		plt.close()
		print("finish drawing accessory ARGs cluster pictures...")

		if len(self.annotation_files) <= 60:
			#draw correlation map
			#set the parameters of the picture
			try:
				plt.rcParams['font.family'] = self.fonttype2
			except:
				print("fonttype not found!")
			plt.rcParams['font.size'] = self.glfs2
			#draw the a blank canvas with number of n^2 subplots
			#(n refers to the number of genomes)
			f,axs = plt.subplots(nrows = subplots_num, ncols = subplots_num,\
				sharex= True, sharey = True, figsize = (self.page_length2,self.page_width2))
			f.subplots_adjust(wspace=0, hspace=0)
			axs[0,0].invert_xaxis()
			#draw the picture in each subplot
			for i in range(subplots_num**2):
				#show the progress of the process
				print("%d/%d drawing No.%d subgraph of total %d subgraphs"\
				 % (i+1,subplots_num**2,i+1,subplots_num**2))
				color_index = []
				#get the position of the subplot
				row_num = i//subplots_num
				col_num = i%subplots_num
				ax = axs[row_num,col_num]
				#give different color of different dot
				for j in range(len(value[row_num])):
					if value[row_num][j] > value[col_num][j]:
						color_index.append("x>y")
					else:
						if value[row_num][j] < value[col_num][j]:
							color_index.append("x<y")
						else:
							color_index.append("x=y")
				#draw the scatterplot map
				sns.scatterplot(x = value[col_num], y = value[row_num], \
					hue = color_index, style = color_index, \
					s= int(self.dotsize), \
					ax = ax, \
					hue_order = ["x=y","x>y","x<y"], style_order = ["x>y","x=y","x<y"], \
					legend = False)
				#rotate the y labels
				if col_num == 0:
					ax.set_ylabel(genome_id[row_num], \
						rotation = 0, rotation_mode='anchor', ha="right", va="top")
				#rotate the x labels
				if row_num == len(genome_id)-1:
					ax.set_xlabel(genome_id[col_num], \
						rotation=40, rotation_mode='anchor', ha="right", va="top")
			#save the correlation map
			print("saving figure...")
			f.savefig(self.analysis_dir+"2_accessory_ar_corr.png", \
				bbox_inches='tight', format="png", dpi=self.dpi2)
			plt.close()
			print("finish drawing accessory ARGs correlation pictures...")
		elif len(self.annotation_files) > 60:
			print("more than 60 genomes are input, correlation picture analysis skipped...")


def main(install_directory,csvfile_directory):
	pansum = Pansum(install_directory,csvfile_directory)
	pansum.pansum_file()
	pansum.pansum_picture()
	if len(pansum.annotation_files) >=2:
		if len(pansum.drug_class) >= 2:
			pansum.pan_picture()
		else:
			print("No more than one drug class was found, pan-resistome analysis skipped...")
	else:
		print("No more than one genome, pan-resistome analysis skipped...")
	accsum = AccessorySummary(pansum.install_dir,pansum.csv_file_dir)
	if len(accsum.annotation_files) >= 2:
		accsum.access_file()
		accsum.accsum_file()
		if len(accsum.drug_title) >= 2:
			accsum.accsum_picture()
			accsum.acc_picture()
		else:
			print("No more than one drug class was found, accessory resistome analysis skipped...")
	else:
		print("No more than one genome, accessory resistome analysis skipped...")
if __name__ == '__main__':
	main("","")