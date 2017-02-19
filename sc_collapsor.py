import os
import sys

def most_common(lst):
    return max(set(lst), key=lst.count)

def determine_tag(file):
	ext = file.split("/")[-1]
	if "v" in ext:
		return "v"
	elif "j" in ext:
		return "j"
	else:
		print "Error: please make sure filename includes v or j for indentification."

def collapse(file):
	print file
	tag = determine_tag(file)
	positions = {} 
	openfile = open(file,"r")
	for line in openfile:
		seq = line.split(", ")[3]
		
		if tag == "j":		# j must be reversed as it aligns from the end of sequences
			seq = seq[::-1]
		for i in range(len(seq)):
			if i in positions:
				positions[i].append(seq[i])
			else:
				positions[i] = []
				positions[i].append(seq[i])
	new_seq_holder = []
	for k in positions:
		new_seq_holder.append(most_common(positions[k]))
	new_seq = ("").join(new_seq_holder)
	if tag == "j":
		new_seq = new_seq[::-1] #j reversed again, to get collapsed seq but in original direction
	print ""
	print new_seq
	print "-------------------------------------------------------------------------------------"
	print ""

directory = sys.argv[1]+"/"
for file in os.listdir(directory):
       collapse(directory+file)