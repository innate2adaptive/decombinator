import os
import io
import sys

def get_files(files_to_split):
	list_of_files=open(files_to_split,"r")
	files = []
	for line in list_of_files:
		files.append(line.rstrip())
	list_of_files.close()
	return files

def tagcollector(filename):
	f = open(filename, 'r')
	v_tags = {}
	j_tags= {}
	double_tags = {}
	for line in f:
		dcr = line.split(", ")
	
		if dcr[0] != "n/a" and dcr[1] == "n/a":
			if dcr[0] in v_tags:
				v_tags[dcr[0]].append(line)
			else:
				v_tags[dcr[0]] = []
				v_tags[dcr[0]].append(line)

		elif dcr[0] == "n/a" and dcr[1] != "n/a":
			if dcr[1] in j_tags:
				j_tags[dcr[1]].append(line)
			else:
				j_tags[dcr[1]] = []
				j_tags[dcr[1]].append(line)

		elif dcr[0] != "n/a" and dcr[1] != "n/a":
			print "--------------------------------------------------------------------------------------------"
			print "double tag found in: "+filename
			print dcr
			print "--------------------------------------------------------------------------------------------"
			if ("v"+dcr[0]+"_j"+dcr[1]) in double_tags:
				double_tags["v"+dcr[0]+"_j"+dcr[1]].append(line)
			else:
				double_tags["v"+dcr[0]+"_j"+dcr[1]] = []
				double_tags["v"+dcr[0]+"_j"+dcr[1]].append(line)

	f.close()
	return v_tags, j_tags, double_tags

def writer(filename, v_tags, j_tags, double_tags):
	directory = filename.split(".")
	directory = "".join(directory[0:len(directory)-1])+"_by_tag"

	if not os.path.exists(directory):
		os.makedirs(directory)

	for v in v_tags.keys(): 
		outfile = directory+"/v"+v+".n12"
		output = open(outfile, "w")
		for dcr_tag in v_tags.get(v):
			output.write(dcr_tag)
		output.close()

	for j in j_tags.keys(): 
		outfile = directory+"/j"+j+".n12"
		output = open(outfile, "w")
		for dcr_tag in j_tags.get(j):
			output.write(dcr_tag)
		output.close()

	for d in double_tags.keys(): 
		outfile = directory+"/"+d+".n12"
		output = open(outfile, "w")
		for dcr_tag in double_tags.get(d):
			output.write(dcr_tag)
		output.close()

files_to_split = sys.argv[1]
files = get_files(files_to_split)
for filename in files:
	tags= tagcollector(filename)
	writer(filename,tags[0],tags[1],tags[2])