
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
# get_barcode_and_UMI7.py
# 1.29.18
#
# Lets try to build this with, you know, functions and things..
# I think this is working. Gives me a vector of dictionary objects, 
# where UMIs are keys and # reads are values. Working - would be 
# better if we could get it to output averages for nReads and nUMIs.   
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
import os
import io
import sys
from shutil import copyfile
import boto3
import botocore
import gzip

#/////////////////////////////////////////////////////////////////////
# dicAppend(): Searches dictionary for item, if found, increases
#               associated value by one, if not, adds to dictionary
#/////////////////////////////////////////////////////////////////////
def dicAppend(d, item):
	if item in d:
		d[item] += 1
	else:
		d[item] = 1
	
#/////////////////////////////////////////////////////////////////////
# find_rUMI(): Searches a string (read) for a substring (Barcode)
#              and then finds the associated right UMI and returns it
#/////////////////////////////////////////////////////////////////////
def find_rUMI(l, s):
	index = l.find(s)
	start = index + 8
	end = index + 53
	u = l[start:end]
	return u 

#/////////////////////////////////////////////////////////////////////
# find_lUMI(): Searches a string (read) for a substring (Barcode)
#               and then finds the associated left UMI, which is 
#               returned
#/////////////////////////////////////////////////////////////////////
def find_lUMI(l, s):
	index = l.find(s)                      
	start = index - 53
	if start < 0:
		start = 0
	end = index - 8
	u = l[start:end]
	return u 

#/////////////////////////////////////////////////////////////////////
# driver(): Main logic of the loop performed here. Give it a file
#           (.fastq) it will search for the associated barcode, 
#           then calls helper functions to add each UMI to a dictionary
#           (cellDic) which is at the end appended to a vector of 
#           dictionary objects (bigList)     
#/////////////////////////////////////////////////////////////////////
def driver(f, b5, b6, b9, b5l, b6l, b9l, nfl, bigList, numReadsList,
                        numUMIList):	

	f1 = f.replace("_R1_001.fastq.gz", "")
        
	#initialize individual cell dictionary
	cellDic = dict()

	b5_found = False
	b6_found = False
	b9_found = False
	
	#loop 1 - find Barcodes
	with io.TextIOWrapper(gzip.open(f, 'r')) as wf:
		all_lines = wf.readlines()
		for line in all_lines:
			if b5 in line:
				b5_found = True
				break
			elif b6 in line:
				b6_found = True
				break
			elif b9 in line:
				b9_found = True
				break
	
	#append to barcodes list
	if b5_found: 
		b5l.append(f1)
	elif b6_found:
		b6l.append(f1)
	elif b9_found:
		b9l.append(f1)
	else:
		nfl.append(f1)

	#loop 2 - find UMIs, add to dictionary
	with io.TextIOWrapper(gzip.open(f, 'r')) as wf:
		all_lines = wf.readlines()
		for line in all_lines:
			if b5 in line:
				u1 = find_lUMI(line, b5)
				u2 = find_rUMI(line, b5)
				dicAppend(cellDic, u1)
				dicAppend(cellDic, u2)					
			elif b6 in line:
				u1 = find_lUMI(line, b6)
				u2 = find_rUMI(line, b6)
				dicAppend(cellDic, u1)
				dicAppend(cellDic, u2)
			elif b9 in line:
				u1 = find_lUMI(line, b9)
				u2 = find_rUMI(line, b9)
				dicAppend(cellDic, u1)
				dicAppend(cellDic, u2)
	wf.close()
	bigList.append(cellDic)
	
	#count reads per cell
	numReads = len(all_lines) / 4
	numReadsList.append(tuple((f1, numReads)))
	
	#count UMIs per cell
	numUMI = len(cellDic)
	numUMIList.append(tuple((f1, numUMI)))
	
	return

#/////////////////////////////////////////////////////////////////////
# getFileNames(): Use boto3 python library to get a list of every file
#					in a specified s3 directory. This list is returned
#					to main function.  
#/////////////////////////////////////////////////////////////////////
def getFileNames(BUCKET_NAME, SUBDIR):

	client = boto3.client('s3')
	paginator = client.get_paginator('list_objects')
	response_iterator = paginator.paginate(Bucket=BUCKET_NAME, 
													Prefix=SUBDIR)

	# this is the list of file names in that bucket starting with that prefix
	file_set = {r['Key'] for result in response_iterator
			for r in result.get('Contents', [])}

	return file_set

#/////////////////////////////////////////////////////////////////////
# downloadFile(): Use boto3 python library to download a single file,
#					the name of which is specified as an arg, from
#					a given s3 bucket/subDir. 
#/////////////////////////////////////////////////////////////////////
def downloadFile(curr_file, file_names_list):
	
	client = boto3.client('s3')
	curr_fileName = curr_file.split('/', 2)[-1]
	
	try:
		client.download_file('lincoln.harris-work', curr_file, curr_fileName)
	except botocore.exceptions.ClientError as e:
		if e.response['Error']['Code'] == "404":
			print("The object does not exist.")
		else:
			raise

	return curr_fileName

#/////////////////////////////////////////////////////////////////////
# main(): Initiates loop that will call driver func for every cell
#          to perform main logic. Initiates lists and globals, opens
#          files, writes to files, etc. Cmd line args:
# 
#           0  - run, but output nothing. Can also be used for testing
#					purposes
#           1  - generate barcode_found output files
# 	    	2  - output # UMIs per well, and average
#           3  - output # reads per well, and average
#           4  - user specifies a cell name, program outputs every UMI
#					for that cell and its associated # of reads
#/////////////////////////////////////////////////////////////////////
	
# initialize some vars
b5="TAGTAGTTCAGACGCCGTTAAGCGC"
b6="CCGTACCTAGATACACTCAATTTGT"
b9="CTGACGTGTGAGGCGCTAGAGCATA"
totalReads = 0
totalUMI = 0

# initialize some lists
b5l = []
b6l = []
b9l = []
nfl = []
bigList = []
numReadsList = []
numUMIList = []

# get a list of every file name in the specified s3 directory
all_file_names = getFileNames('lincoln.harris-work','171128_M05295_0059_000000000-BFG9N')

# filter, so that we're left with just the fastq files in that s3 dir
fastq_file_names = []
for item in all_file_names:
	if('.fastq.gz' in item):
		fastq_file_names.append(item)

# MAIN LOOP 
# 	For each iteration of loop: 
#		download (single) .fastq file
#		call driver func to create cell dictionary, and get UMI and read counts
for item in fastq_file_names:
	f_name = downloadFile(item, fastq_file_names)
	driver(f_name, b5, b6, b9, b5l, b6l, b9l, nfl, bigList, numReadsList, 
			numUMIList)
	os.remove(f_name)

# write lists to files    works!!!
if sys.argv[1] == "1":
	b5f = open("b5_found.txt", "w")
	b6f = open("b6_found.txt", "w")
	b9f = open("b9_found.txt", "w")
	nf = open("notFound.txt", "w")
	for item in b5l:
		b5f.write("%s\n" % item)
	for item in b6l:
		b6f.write("%s\n" % item)
	for item in b9l:
		b9f.write("%s\n" % item)
	for item in nfl:
		nf.write("%s\n" % item)
	b5f.close()
	b6f.close()
	b9f.close()
	nf.close()

for item in numUMIList:
	totalUMI += item[1]

for item in numReadsList:
	totalReads += item[1]


# find num UMIs per well
if sys.argv[1] == "2":
	with open("numUMIs_out.txt", "w") as of:
		of.write('\n'.join('%s %d' % x for x in numUMIList))

	aveUMI = totalUMI/len(numUMIList)
	print(" ")
	print("####################################################")
	print("         Average number of UMIs per well: ")
	print("                   %d " % aveUMI)
	print("####################################################")
	print(" ")

# find num reads per well
if sys.argv[1] == "3":
	with open("numReads_out.txt", "w") as of:
		of.write('\n'.join('%s %d' % x for x in numReadsList))		
		
	aveReads = totalReads/len(bigList)
	print(" ")
	print("####################################################")
	print("         Average number of reads per well: ")
	print("                   %d " % aveReads)
	print("####################################################")
	print(" ")

# display UMIs and number of associated reads for each one, for a specified cell
if sys.argv[1] == "4":
	cellList = []
	for item in fastq_file_names:
		item1 = item.split('/', 2)[-1]
		item2 = item1.replace("_R1_001.fastq.gz", "")
		cellList.append(item2)

	if len(sys.argv) != 3:
		print(" ")
		print("ERROR. Enter a valid cell name.")
		print(" ")
		sys.exit()
	
	qCell = sys.argv[2]
	if qCell not in cellList:
		print(" ")
		print("ERROR. Enter a valid cell name.")
		print(" ")
		sys.exit()

	index = cellList.index(qCell)
	
	uList = []
	
	#send to an out file
	for item in bigList[index]:
		tup = (item, bigList[index][item])
		uList.append(tup)
		#print(item, bigList[index][item])  #works like a charm!!
	
	with open("UMIs_out.txt", "w") as of:
                of.write('\n'.join('%s %d' % x for x in uList))

	print(" ")
	print("cell ID: %s" % qCell)
	print("total reads: %d" % totalReads)
	print("total UMIs: %d" % totalUMI)
	print(" ")


#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////
#/////////////////////////////////////////////////////////////////////

