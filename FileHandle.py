import os,csv,sys
#The moudle is written to handle the files and directories

#get all the files with a given extension in a given directory
def filename_get(directory,extension,showext=False):
	file_name = []
	filenames = os.listdir(directory)
	for filename in filenames:
		if filename.endswith(extension):
			if not showext:
				file_name.append(filename.split(str(extension))[0])
			else:
				file_name.append(filename)
	return file_name

#get the current working directory
def cwd_get(directory):
	if directory == "":
		file_directory = (os.getcwd()).replace("\\","/")+"/"
		#use "/" to separate the path
	else:
		file_directory = directory
	return file_directory

#get the input directory
def inpd_get(directory="",file_categories=""):
	#if no directory is input, then ask for the required directory
	if directory == "":
		file_directory = (input("please input the %s directory: " % file_categories)).replace("\\","/")+"/"
	#if directory is input, return the input directory
	else:
		file_directory = directory
	return file_directory

#create a new directory if the directory not exists
def dir_add(file_directory):
	if not os.path.exists(file_directory):
		os.mkdir(file_directory)

#remove a file if it exists
def file_del(filename):
	if os.path.exists(filename):
		os.remove(filename)

#add each line of the file into a list
def file_reader(filename):
	text = []
	with open(filename,"r") as f:
		text = list(f.readlines())
	return text

#read the csv extention files
def csv_reader(csvfilename,delim=","):
	plots = []
	with open(csvfilename,"r") as csvfile:
		plots = list(csv.reader(csvfile, delimiter=delim))
		#retun the file content in a list
		#each sublist in the list represent a row in the sheet
	return plots

def setting_reader(install_directory,setting_name):
	set_return = ""
	with open(install_directory+"settings.txt","r") as set_file:
		settings = list(set_file.readlines())
		for line in settings:
			if setting_name in line and "[" in line:
				set_return = line.split("=")[1].split("]")[0]
				if "blast" in setting_name:
					if not set_return.endswith("/") or set_return.endswith("\\"):
						set_return = set_return + "/"
	return set_return

#write the content into a file
def file_write(filename,content):
	with open(filename,"a") as f:
		f.write(content)
