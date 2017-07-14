import sys

def get_files(input_file,output_file):
	files_to_merge = []
	merged_file_names = []
	
	open_input_file = open(input_file,"r")
	for line in open_input_file:
		files_to_merge.append(line.rstrip())
	open_input_file.close()

	open_output_file = open(output_file,"r")
	for line in open_output_file:
		merged_file_names.append(line.rstrip())
	open_output_file.close()

	return files_to_merge, merged_file_names


def writer(files_to_merge,output_files):
	output_count = 0	

	for i in range(0,8,2):
		print "--------------------------"
		print ""
		print "files to merge:"
		print ""
		print files_to_merge[i]
		print files_to_merge[i+1]
		print ""
		
		file1 = open(files_to_merge[i],"r")
		file1_lines = []
		for line in file1:
			file1_lines.append(line)
		file1.close()
		
		file2 = open(files_to_merge[i+1],"r")
		file2_lines = []
		for line in file2:
			file2_lines.append(line)
		file2.close()	

		outfile = open(output_files[output_count],"w")
		for i in file1_lines:
			outfile.write(i)
		for i in file2_lines:
			outfile.write(i)
		outfile.close()	

		print "data written to ", output_files[output_count]
		print ""
		output_count+=1

input_file_names = sys.argv[1]
output_file_names = sys.argv[2]

file_names = get_files(input_file_names,output_file_names)
writer(file_names[0],file_names[1])