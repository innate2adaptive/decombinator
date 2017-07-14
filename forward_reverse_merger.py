import sys

def gather_dcrs(a_direction, b_direction):

	a_dcrs =[]
	a_file = open(a_direction, 'r')
	for line in a_file:
		a_dcrs.append(line)
	a_file.close()	

	b_dcrs =[]
	b_file = open(b_direction, 'r')
	for line in b_file:
		b_dcrs.append(line)
	b_file.close()

	'''
	-------------------------------------------------------------------------------------------------
	-------------------------------------------------------------------------------------------------
	The following loop was a test to find any duplicate results in reverse and forward directions.
	No duplicates were found in any paired files, so the files can be simply merged. If duplicates were
	found, a small amount of additional functionality would need to be added to the writer function, 
	to prevent these duplicates being written to the output files.
	-------------------------------------------------------------------------------------------------

	for a in a_dcrs:
		a_fastqid = a.split(',')[2]
		a_seq = a.split(',')[3]
		a_qual = a.split(',')[4]
		for b in b_dcrs:
			b_fastqid = b.split(',')[2]
			b_seq = b.split(',')[3]
			b_qual = b.split(',')[4]
			if (b_fastqid == a_fastqid):
				print "-----------------------------------------------------------------"
				print "match!"
				print "a: "+a
				print "b: "+b
				print "-----------------------------------------------------------------"

-	------------------------------------------------------------------------------------------------
	-------------------------------------------------------------------------------------------------
	'''
	return a_dcrs, b_dcrs

def writer(a_dcrs, b_dcrs, outfile):
	output = open(outfile, "w")
	print "Files merged to " + outfile
	for a in a_dcrs:
		output.write(a)
	for b in b_dcrs:
		output.write(b)
	output.close()

def file_namer(a_file, b_file):
	a_file = a_file.split('_')
	b_file = b_file.split('_')
	del a_file[-1]
	del b_file[-1]
	if a_file != b_file:
		print "Error: files from different sources. Make sure paired files have matching read number and chain."
	else:
		prefix = '_'.join(a_file)
		merge_name = prefix + "_merge.n12"

		return merge_name

def produce_merged_files(list_of_files):
	files_to_merge = open(list_of_files,"r")
	files = []
	for line in files_to_merge:
		files.append(line.rstrip())
	files_to_merge.close()

	for i in range(0,len(files),2):
		forward = files[i]
		reverse = files[i+1]

		merged_file_name = file_namer(forward,reverse)
	
		dcrs = gather_dcrs(forward,reverse)
		writer(dcrs[0],dcrs[1],merged_file_name)


list_of_files = sys.argv[1]
produce_merged_files(list_of_files)