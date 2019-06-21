import os,csv,sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns
import random
import itertools
import FileHandle as fh
from scipy.optimize import curve_fit
from scipy.special import comb
from decimal import *

class Distribution:

	#attributes of class Distribution
	#the directory of software, csv annotation files and pangenome analysis files
	install_dir = ""
	csv_file_dir = ""
	pangenome_dir = ""
	#parameters for drawing picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi = 200
	fontsize = 20
	fonttype = "Times New Roman"
	page_length = 15
	page_width = 15
	annotation_files = []

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
		self.pangenome_dir = self.csv_file_dir+"pangenome/"
		fh.dir_add(self.pangenome_dir)
		#get the parameters for the picture from settings.txt (1)ar distribution picture
		self.dpi = int(fh.setting_reader(self.install_dir,"dpi1"))
		self.fontsize = int(fh.setting_reader(self.install_dir,"fontsize1"))
		self.fonttype = fh.setting_reader(self.install_dir,"fonttype1")
		self.page_length = int(fh.setting_reader(self.install_dir,"page_length1"))
		self.page_width = int(fh.setting_reader(self.install_dir,"page_width1"))

	#count all the ar genes in the annotation files and return a dictionary
	def ar_count(self):
		#a list to store all the ar genes
		pan_ar_gene_sum = []
		#a dict which keys represent ar gene names and values represent times appeared
		pan_ar_gene_count = {}
		for each in self.annotation_files:
			#store the annotation information into a list
			plots = fh.csv_reader(self.csv_file_dir+each)
			for row in plots[1:]:
				pan_ar_gene_sum.append([row[0],row[6],row[2]])
		#count the ar genes and write the result into a dictionary
		for each in pan_ar_gene_sum:
			if each[1] not in pan_ar_gene_count:
				pan_ar_gene_count[each[1]] = 1
			else:
				pan_ar_gene_count[each[1]] += 1
		return pan_ar_gene_count

	#use the dictionary from ar_count() to write a file and return a dataframe to draw picture
	def distribution_file(self,ar_gene_dict):
		print("ar distribution analysis start")
		#create the file under the pangenome result directory
		ar_distribution_file = self.pangenome_dir + "4_ar_distribution.txt"
		fh.file_del(ar_distribution_file)
		f_ar_distribution = open(ar_distribution_file,"a")
		ar_distribution = {}
		#classify the ar genes according to there distribution number
		for count in range(1,len(self.annotation_files)+1):
			line = str(count)+"	"
			ar_distribution[count] = []
			for key in ar_gene_dict:
				if ar_gene_dict[key] == count:
					ar_distribution[count].append(key)
					line += key+","
			line = line.strip(",")+"\n"
			f_ar_distribution.write(line)
		f_ar_distribution.close()
		return ar_distribution

	#draw the distribution picture
	def distribution_picture(self,ar_distribution_dict):
		#determine the fontsize and fonttype according to settings
		try:
			plt.rcParams['font.family'] = self.fonttype
		except:
			print("fonttype not found!")
		plt.rcParams['font.size'] = self.fontsize
		#x labels are number of genomes
		ar_distri_xlable = [key for key in ar_distribution_dict.keys()]
		#values are number of genes in the certain number of genomes
		ar_distri_value = [len(ar_distribution_dict[key]) for key in ar_distribution_dict.keys()]
		#determine the figure size according to settings
		f,ax = plt.subplots(figsize = (self.page_length,self.page_width))
		#draw the barplot of ar genes distribution
		sns.barplot(x = ar_distri_xlable, y = ar_distri_value, palette = "Blues_d")
		ax.set_xticks([0,int(len(self.annotation_files)/5)-1,int(len(self.annotation_files)/5)*2-1,\
			int(len(self.annotation_files)/5)*3-1,int(len(self.annotation_files)/5)*4-1,int(len(self.annotation_files)/5)*5-1])
		ax.set_xticklabels([1,int(len(self.annotation_files)/5),int(len(self.annotation_files)/5)*2,\
			int(len(self.annotation_files)/5)*3,int(len(self.annotation_files)/5)*4,int(len(self.annotation_files)/5)*5])
		#set name of x axis
		ax.set_xlabel('Number of Genomes')
		#set name of y axis
		ax.set_ylabel('Number of Antibiotic Resistance Genes')
		#save figure in png format and dpi in settings
		f.savefig(self.pangenome_dir+"4_ar_distribution.png", bbox_inches='tight', format="png", dpi=self.dpi)
		plt.close()
		print("ar distribution analysis finished...")

class Pangenome:
	#attributes of class Distribution
	#the directory of software, csv annotation files and pangenome analysis files
	install_dir = ""
	csv_file_dir = ""
	pangenome_dir = ""
	#parameters for drawing picture
	#dpi(dot per inch),font size and type, picture length and width
	dpi = 200
	fontsize = 20
	fonttype = "Times New Roman"
	page_length = 20
	page_width = 15
	annotation_files = []
	#the parameters of fitting model and coverage
	fit_coverage = 0.8
	fit_model = "power_law"

	def __init__(self,install_directory,csvfile_directory):
		#get the directory of software, csv annotation files
		self.install_dir = fh.cwd_get(install_directory)
		self.csv_file_dir = fh.inpd_get(csvfile_directory,file_categories="csv files")
		#create the directory of pangenome analysis files
		self.pangenome_dir = self.csv_file_dir+"pangenome/"
		fh.dir_add(self.pangenome_dir)
		self.pangenome_file = self.pangenome_dir + "5_pangenome.txt"
		self.fit_model_file = self.pangenome_dir+"5_pangenome_fitmodel.txt"
		#add all annotation files name into a list with extension
		self.annotation_files = fh.filename_get(self.csv_file_dir,"_ar.csv",showext=True)
		#get the parameters for the picture from settings.txt (1)ar distribution picture
		self.dpi = int(fh.setting_reader(self.install_dir,"dpi2"))
		self.fontsize = int(fh.setting_reader(self.install_dir,"fontsize2"))
		self.fonttype = fh.setting_reader(self.install_dir,"fonttype2")
		self.page_length = int(fh.setting_reader(self.install_dir,"page_length2"))
		self.page_width = int(fh.setting_reader(self.install_dir,"page_width2"))
		#get the parameters of fitting model and coverage
		self.fit_model = fh.setting_reader(self.install_dir,"fit_model")
		self.fit_coverage = float(fh.setting_reader(self.install_dir,"fit_coverage"))

	def pan_file(self):

		#calculate combinations in C(k,N), return a list of combinations.
		#number of the combination in combinations is equal to max_combs
		def Combinations(iterable,num_elements,max_combs=10000):
			all_combinations = []
			#ensure the same random combination formed
			random.seed(0)
			#only when combination exceed twice of the max_combs
			#the combination will be decreased to max_combs
			#this is to save time while exclude the same random combination 
			#based on this algorithm
			if comb(len(iterable),num_elements) > max_combs*2:
				#generate different combinations until the number reach max_combs
				while len(all_combinations) < max_combs:
					#generate a random combination from the iterable list
					each_combination = random.sample(iterable,num_elements)
					#add the different combination into the result list
					if each_combination not in all_combinations:
						all_combinations.append(each_combination)
			else:
				#if number of combination is smaller than twice max_combs
				#return all the possible combinations
				all_combinations = list(itertools.combinations(iterable,num_elements))
				#result may also exceed max_combs
				#so just randomly select max_combs of them from the result
				if len(all_combinations) > max_combs:
					all_combinations = random.sample(all_combinations,max_combs)
			#return all the combinations
			return all_combinations

		print("pangenome analysis start...")
		#store all the annotation information into a dictionaty to save time on reading
		filename_arlist = {}
		for each in self.annotation_files:
			filename_arlist[each] = fh.csv_reader(self.csv_file_dir+each)
		#choose to run pangenome analysis by default
		choose = "Y"
		#if pangenome file already exists, you can choose whether to update or not
		#this step may waste time, you can skip if you just want to re-fit the model
		#without change of original data
		if os.path.exists(self.pangenome_file):
			choose = input("update pangenome file (Y/N)?:").upper()
		if choose == "Y":
			fh.file_del(self.pangenome_file)
			f_pangenome = open(self.pangenome_file,"a")
			#the title of the pangenome file
			#the first row is genome number, the second is pangenome number
			#and the last is core genome number
			f_pangenome.write("Genome_Number"+"	"+"Pan_Genome_Size"+"	"+"Core_Genome_Size\n")
			#calculate the pan and core genome size under each number of combination
			for i in range(1,len(self.annotation_files)+1):
				print("calculating combiantions in C("+str(len(self.annotation_files))+","+str(i)+")...")
				#get the list of possible combinations
				combinations = Combinations(self.annotation_files,i)
				for combination in combinations:
					#dictionary to store all ar genes and there amounts
					comb_ar_gene_all = {}
					#list to store ar genes in pangenome
					comb_ar_gene_pan = []
					#list to store ar genes in core genome
					comb_ar_gene_core = []
					#traverse all the genomes in the combination
					for file in combination:
						file_content = filename_arlist[file]
						#count the amount of ar genes
						for row in file_content:
							if row[1] != "gene_name":
								if row[6] not in comb_ar_gene_all.keys():
									comb_ar_gene_all[row[6]] = 1
								else:
									comb_ar_gene_all[row[6]] += 1
					#classify the ar genes into pan and core according to their amounts
					for gene in comb_ar_gene_all.keys():
						comb_ar_gene_pan.append(gene)
						if comb_ar_gene_all[gene] >= len(combination):
							comb_ar_gene_core.append(gene)
					#write the amounts of genomes, pan ar genes and core ar genes into the file
					f_pangenome.write(str(len(combination))+"	"+str(len(comb_ar_gene_pan))+"	"+str(len(comb_ar_gene_core))+"\n")
			f_pangenome.close()
		else:
			print("calculating pangenome skipped...")
		print("pangenome analysis finished...")
	
	def pan_fit(self):
		print("draw picture start...")
		#determine the fontsize and fonttype according to settings
		try:
			plt.rcParams['font.family'] = self.fonttype
		except:
			print("fonttype not found!")
		plt.rcParams['font.size'] = self.fontsize
		#get x value and y1,y2 values from pangenome file
		plots = fh.csv_reader(self.pangenome_file, delim="	")
		x_value = [int(row[0]) for row in plots[1:]]
		y1_value = [int(row[1]) for row in plots[1:]]
		#define the category of y1 value
		y1_categories = ["Pan Resistome" for row in plots[1:]]
		y2_value = [int(row[2]) for row in plots[1:]]
		#define the category of y2 value
		y2_categories = ["Core Resistome" for row in plots[1:]]
		#create the dataframe of the pangenome
		#including x-y vaules and thier lables related to their categories
		data = {"Number of genomes":x_value+x_value,"Number of Antibiotic Resistance Genes":y1_value+y2_value,"Genome_categories":y1_categories+y2_categories}
		df = pd.DataFrame(data)
		f,ax = plt.subplots(figsize = (self.page_length,self.page_width))
		#draw boxplot pictures to show the distribution of each genome size on the picture
		sns.boxplot(data = df, x = "Number of genomes", y = "Number of Antibiotic Resistance Genes", hue = "Genome_categories", width = 0.5, linewidth = 1, fliersize = 1, dodge = False, ax = ax)

		#average values of genome amounts
		x_average = []
		#average values of x-1
		x_average_minus = []
		#average values of pangenome size in each amount of genome
		y1_average = []
		#average values of core genome size in each amount of genome
		y2_average = []
		y_index_start = 0
		y_index_end = 0
		#calculate the average values according to the data
		for i in range(len(x_value)):
			#get all the amount of genomes
			if x_value[i] not in x_average:
				x_average.append(x_value[i])
				x_average_minus.append(x_value[i]-1)
				y_index_end = i
				#sum all the pan and core genome value in a fixed amount if genome
				if y_index_start < y_index_end:
					y1_value_sum = 0
					y2_value_sum = 0
					for m in range(y_index_start,y_index_end):
						y1_value_sum += y1_value[m]
						y2_value_sum += y2_value[m]
					y1_average.append(y1_value_sum/(y_index_end-y_index_start))
					y2_average.append(y2_value_sum/(y_index_end-y_index_start))
					y_index_start = y_index_end
			#the last amount of the genome may be ignored by the method
			#so the last value of pan and core genome is append to the list
			if i == len(x_value)-1:
				y1_average.append(y1_value[i])
				y2_average.append(y2_value[i])

		#use numpy.array(list) to transfer the list into array format
		x_value = np.array(x_value)
		y1_value = np.array(y1_value)
		y2_value = np.array(y2_value)
		x_average = np.array(x_average)
		x_average_minus = np.array(x_average_minus)
		y1_average = np.array(y1_average)
		y2_average = np.array(y2_average)

		#draw the curve of average values in the format of "-----"
		ax.plot(x_average_minus, y1_average, color='lightblue', lw=2.0, ls="--", label = "Curve of Average Pan Resistome Size")
		ax.plot(x_average_minus, y2_average, color='orange', lw=2.0, ls="--", label = "Curve of Average Core Resistome Size")

		#the percentage that will not be calculated in the fitting
		notcal = 1-self.fit_coverage
		genome_num = len(x_average)
		#cut the list according to the fitting coverage
		cut = int(genome_num*notcal)
		x_average = x_average[cut:]
		x_average_minus = x_average_minus[cut:]
		y1_average = y1_average[cut:]
		y2_average = y2_average[cut:]

		#used to calculate the R square according to the fitting model
		def Rcal(y_average,y_fit):
			genome_num = len(x_average)
			y_sum = 0
			for i in range(genome_num):
				y_sum += y_average[i]
			y_mean = y_sum/genome_num
			SS1 = Decimal(0)
			SS2 = Decimal(0)
			for i in range(genome_num):
				SS1 += (Decimal(y_average[i])-Decimal(y_fit[i]))**2
				SS2 += (Decimal(y_average[i])-Decimal(y_mean))**2
			R_square = 1-Decimal(SS1)/Decimal(SS2)
			return R_square

		#create the file to save fitting models
		fh.file_del(self.fit_model_file)
		f_fit_model = open(self.fit_model_file,"a")

		#if the fitting model is "polyfit" in settings.txt
		if self.fit_model == "polyfit":

			#used to write the model into the file according to different fitting order
			def polyfit_write(fit_list,fit_order):
				for i in range(fit_order+1):
					#while x in each item
					if i < fit_order:
						f_fit_model.write("("+str(fit_list[i])+")*x^"+str(fit_order-i)+"+")
					#the constant term
					elif i == fit_order:
						f_fit_model.write("("+str(fit_list[i])+")")
			#get the fitting order from settings.txt
			ployfit_order = int(fh.setting_reader(self.install_dir,"fit_order"))
			#use numpy.polyfit to fit the model
			poly1 = np.polyfit(x_average, y1_average, ployfit_order)
			#calculate the fit y value according to the model
			y1_fit= np.polyval(poly1, x_average)
			#draw the fitting curve of the model
			ax.plot(x_average_minus, y1_fit, color='purple', lw=3.0, label = "Fitting Curve of Pan Resistome Size")
			
			poly2 = np.polyfit(x_average, y2_average, ployfit_order)
			y2_fit= np.polyval(poly2, x_average)
			ax.plot(x_average_minus, y2_fit, color='red', lw=3.0, label = "Fitting Curve of Core Resistome Size")
			#calculate the R square
			R_square1 = Rcal(y1_average,y1_fit)
			R_square2 = Rcal(y2_average,y2_fit)
			#write down the fitting model
			f_fit_model.write("Polynomial Regression Model of Pan Resistome:\nP=")
			polyfit_write(poly1,ployfit_order)
			f_fit_model.write(" (R^2="+str(R_square1)+")\n")

			f_fit_model.write("Polynomial Regression Model of Core Resistome:\nC=")
			polyfit_write(poly2,ployfit_order)
			f_fit_model.write(" (R^2="+str(R_square2)+")\n")
			f_fit_model.close()

		#if the fitting model is "power_law" in settings.txt
		elif self.fit_model == "power_law":
			#define the power_law modle
			def power(x,a,b):
				return a*pow(x,b)
			#there may be no proper parameters, so RuntimeError will be raised
			try:
				#power_law model of the pangenome
				fit_res1,pocv = curve_fit(power,x_average,y1_average)
				y1_fit = [power(x,fit_res1[0],fit_res1[1]) for x in x_average]
				ax.plot(x_average_minus, y1_fit, color='purple', lw=3.0, label = "Fitting Curve of Pan Resistome Size")
				#power_law fit of the core genome
				fit_res2,pocv = curve_fit(power,x_average,y2_average)
				y2_fit = [power(x,fit_res2[0],fit_res2[1]) for x in x_average]
				ax.plot(x_average_minus, y2_fit, color='red', lw=3.0, label = "Fitting Curve of Core Resistome Size")
				#calculate the R square
				R_square1 = Rcal(y1_average,y1_fit)
				R_square2 = Rcal(y2_average,y2_fit)
				#write down the fitting model
				f_fit_model.write("Power Law Model of Pan Resistome:\nP=("\
					+str(fit_res1[0])+")*x^("+str(fit_res1[1])+") (R^2="+str(R_square1)+")\n")
				f_fit_model.write("Power Law Model of Core Resistome:\nC=("\
					+str(fit_res2[0])+")*x^("+str(fit_res2[1])+") (R^2="+str(R_square2)+")\n")
			#raise RuntimeError
			except RuntimeError:
				f_fit_model.write("No model fitting")
				print("Optimal parameters not found!")

		elif self.fit_model == "pangp":
			#define the power_law modle
			def power(x,a,b,c):
				return a*pow(x,b)+c
			def exp(x,a,b,c):
				return a*np.exp(x*b)+c
			#there may be no proper parameters, so RuntimeError will be raised
			try:
				#power_law model of the pangenome
				fit_res1,pocv = curve_fit(power,x_average,y1_average)
				y1_fit = [power(x,fit_res1[0],fit_res1[1],fit_res1[2]) for x in x_average]
				ax.plot(x_average_minus, y1_fit, color='purple', lw=3.0, label = "Fitting Curve of Pan Resistome Size")
				#power_law fit of the core genome
				fit_res2,pocv = curve_fit(power,x_average,y2_average)
				y2_fit = [power(x,fit_res2[0],fit_res2[1],fit_res2[2]) for x in x_average]
				ax.plot(x_average_minus, y2_fit, color='red', lw=3.0, label = "Fitting Curve of Core Resistome Size")
				#calculate the R square
				R_square1 = Rcal(y1_average,y1_fit)
				R_square2 = Rcal(y2_average,y2_fit)
				#write down the fitting model
				f_fit_model.write("Power Law Model of Pan Resistome:\nP=("\
					+str(fit_res1[0])+")*x^("+str(fit_res1[1])+")+("+str(fit_res2)+") (R^2="+str(R_square1)+")\n")
				f_fit_model.write("Power Law Model of Core Resistome:\nC=("\
					+str(fit_res2[0])+")*e^("+str(fit_res2[1])+"*x)+("+str(fit_res2)+") (R^2="+str(R_square2)+")\n")
			#raise RuntimeError
			except RuntimeError:
				f_fit_model.write("No model fitting")
				print("Optimal parameters not found!")
		elif self.fit_model == "False":
			f_fit_model.write("No model fitting")

		else:
			print("please input the right model name!")
		#put the lengend to center right
		plt.legend(loc="center right")
		ax.set_xticks([0,int(len(self.annotation_files)/5)-1,int(len(self.annotation_files)/5)*2-1,\
			int(len(self.annotation_files)/5)*3-1,int(len(self.annotation_files)/5)*4-1,int(len(self.annotation_files)/5)*5-1])
		ax.set_xticklabels([1,int(len(self.annotation_files)/5),int(len(self.annotation_files)/5)*2,\
			int(len(self.annotation_files)/5)*3,int(len(self.annotation_files)/5)*4,int(len(self.annotation_files)/5)*5])
		#save the figure
		f.savefig(self.pangenome_dir+"5_pangenome.png", bbox_inches='tight', format="png", dpi=self.dpi)
		plt.close()

		print("draw picture finished...")
		print("pangenome analysis finished...")

def main(install_directory,csvfile_directory):
	dis = Distribution(install_directory,csvfile_directory)
	ar_gene_dict = dis.ar_count()
	ar_distribution_dict = dis.distribution_file(ar_gene_dict)
	dis.distribution_picture(ar_distribution_dict)
	pan = Pangenome(dis.install_dir,dis.csv_file_dir)
	pan.pan_file()
	pan.pan_fit()

if __name__ == '__main__':
	main("","")