#! /usr/bin/python3

# Written by D. Brink 2015-11-22. Last updated 2016-05-04.
# For use with .fcs files exported from the BD Accuri C6 software
#
# INPUT:text file with names of each well in the format ROW_<POS>:<name> e.g. 
#		'C05:glucose 20g.l rep1' and each column as COL_<column digit>:<name> e.g.COL_01:strain1 rep1
#		Input file must also include paths to the .fcs files that should be renamed. 
#		Optionally, renamed files can be moved to new paths, and directories can be created, by stating such in the input text file


import sys, re, os, shutil


## Check if all input files are given ##
if len(sys.argv)>1:
	f_in = open(sys.argv[1], 'r')
else:
	sys.exit('ERROR: Please input plate matrix')

print('Please Wait.')

## Read the input ##
plate_names=[]
plate_dict={}
rename_count=0
maindir=[]
subdir=[]
sub_subdir=[]
row_dict={}

for line in f_in:
	if line !='\n' and not line.startswith('#'):	#excludes empty lines (e.g. final line)
		line=line.rstrip()
		if line.startswith('PATH='): 
			file_dir=line.replace('PATH=','')
		elif line.startswith('NEWPATH='):		
			new_file_dir=line.replace('NEWPATH=','')
		elif line.startswith('MAIN_DIR='):
			line=line.replace('MAIN_DIR=','')
			maindir.append(line)
		elif line.startswith('SUBDIR='):
			line=line.replace('SUBDIR=','')
			subdir.append(line)						#must use append here to put the whole string in an element in the list
		elif line.startswith('SUB_SUBDIR='):		
			line=line.replace('SUB_SUBDIR=','')
			sub_subdir.append(line)
		elif line.startswith('COL_'):
			line=line.replace('COL_','')
			for row_size in range(0, len(plate_names)):
					row_name_and_string=re.split(':',line)
					if len(row_name_and_string)==2:
						row_dict[row_name_and_string[0]]=row_name_and_string[1]
		elif line.startswith('ROW_'):
			line=line.replace('ROW_','')
			plate_names=re.split('\t',line)			# if user has entered the names in a tsv format instead of line by line
			for row_size in range(0, len(plate_names)):
				name_and_string=re.split(':',plate_names[row_size])
				if len(name_and_string)==2:
					plate_dict[name_and_string[0]]=name_and_string[1]
				else:
					sys.exit('ERROR: Each entry in input plate matrix must start with ROW_POS:, e.g. A01: gluc 20 g.l rep1')


## Copy the .fcs files to a temp directory and rename them according to the plate_names list ##
if file_dir:
	os.chdir(file_dir)
else:
	sys.exit('ERROR:No path to the files that should be renamed was given in the input file')
if not os.path.exists('./raw_data'):
	os.mkdir('raw_data')

if new_file_dir !='':
	if not os.path.exists(new_file_dir):
		os.mkdir(new_file_dir)
elif new_file_dir =='':
	new_file_dir=file_dir

source=os.listdir('.')
destination='./raw_data/'
for files in source:
	if files.endswith('.fcs'):
		shutil.copy(files, destination)

row_names=['A','B','C','D','E','F','G','H']
column_names=['01','02','03','04','05','06','07','08','09','10','11','12']
time_dir=[]

for i_rows in range(0, len(row_names)):
	for i_columns in range(0,len(column_names)):
		plate_well_name=row_names[i_rows]+column_names[i_columns]
		for filename in os.listdir('.'):
			curr_dir=os.getcwd()
			if re.search('\d+h',curr_dir, re.I):
				time_dir=re.search('\d+h',curr_dir, re.I).group(0) 	#.group(0) extracts the match to var
			if column_names[i_columns] in row_dict:
				extra_info= ' ' + row_dict[column_names[i_columns]] + ' ' + time_dir
			if filename.startswith(plate_well_name) and filename.endswith('.fcs') and plate_well_name in plate_dict:
				os.rename(filename, plate_dict[plate_well_name] + extra_info + '.fcs')	
				rename_count+=1
				if new_file_dir is not file_dir:		
					shutil.move(plate_dict[plate_well_name] + extra_info + '.fcs', new_file_dir)

print('The script has renamed', rename_count, '.fcs files. Now moving files to desired directories.')


### Sort renamed files for use with the FI_size_normalization.m script  
if new_file_dir:
	os.chdir(new_file_dir)
else:
	os.chdir(file_dir)

# Move files to the designated main directories
if maindir:
	for renamed_file in os.listdir('.'):
		renamed_file=renamed_file.rstrip()
		for j_maindir in range(0,len(maindir)):
			if re.search(maindir[j_maindir],renamed_file, re.I):
				if not os.path.exists(new_file_dir + '/' + maindir[j_maindir]):
					os.mkdir(new_file_dir + '/' + maindir[j_maindir])
				shutil.move(renamed_file,new_file_dir + '/' + maindir[j_maindir])

# Move files to the designated subdirectories
if maindir and subdir:
	for i_maindir in range(0,len(maindir)):
		os.chdir(new_file_dir + '/' + maindir[i_maindir])		#in turn, go to the different subdirs that were created by the last code block
		for renamed_file in os.listdir('.'):
			renamed_file=renamed_file.rstrip()
			for j_subdir in range(0,len(subdir)):
				if re.search(subdir[j_subdir],renamed_file, re.I):
					if not os.path.exists('./' + subdir[j_subdir]):
						os.mkdir(subdir[j_subdir])
					shutil.move(renamed_file,'./' + subdir[j_subdir])

# Move files to the designated sub subdirectories
if maindir and subdir and sub_subdir:
	for i_maindir in range(0,len(maindir)):
		os.chdir(new_file_dir + '/' + maindir[i_maindir])
		for k_subdir in range(0,len(subdir)):
			os.chdir(new_file_dir + '/' + maindir[i_maindir] +'/' + subdir[k_subdir])	#in turn, go to the different sub-subdirs that were created by the last code block
			for subdir_files in os.listdir('.'):
				subdir_files=subdir_files.rstrip()
				for j_sub_subdir in range(0,len(sub_subdir)):
					if re.search('\s'+ sub_subdir[j_sub_subdir]+'\s',subdir_files,re.I):#if re.search(sub_subdir[j_sub_subdir]+'\s',subdir_files,re.I):		#OLD CODE: if re.search('^' + sub_subdir[j_sub_subdir]+'\s',subdir_files,re.I):
						if not os.path.exists('./' + sub_subdir[j_sub_subdir]):
							os.mkdir(sub_subdir[j_sub_subdir])
						shutil.move(subdir_files,'./' + sub_subdir[j_sub_subdir])

print('Run completed. Thank you.')

