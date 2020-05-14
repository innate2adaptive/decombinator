import os
import sys
import subprocess
import argparse

def args():
	parser = argparse.ArgumentParser(
		description='Run pipeline on the test data')
	parser.add_argument('-c1', '--checkpoint1', type=str, required=False,
		help='Choose which script is first to run out of pipeline',
		default="demultiplexor")
	parser.add_argument('-c2', '--checkpoint2', type=str, required=False,
		help='Choose which script is last to run out of pipeline',
		default="cdr3translator")
	return parser.parse_args()	

def getfiles(search, filedir):
	listdir = os.listdir(filedir)
	
	if type(search) == str:
		search = [search]

	foundfiles = []
	for fname in listdir:
		found = True
		for i in search:
			if i.lower() not in fname.lower():
				found = False
		if found:
			foundfiles.append(filedir+os.sep+fname)
	# results = filter(lambda fname: search in fname.lower(), listdir)

	return foundfiles

def execute(command):
	process = subprocess.Popen(command, stdout=subprocess.PIPE)
	output, error = process.communicate()
	print(output.decode())
	return output, error

def demultiplex(test_data_dir):
	## Demultiplexor
	r1 = getfiles("r1", test_data_dir)
	r2 = getfiles("r2", test_data_dir)
	i1 = getfiles("i1", test_data_dir)
	ind = getfiles("indexes", test_data_dir)
	demultfiles = {"r1": None,"r2": None, "i1": None,"ind": None}

	for search in demultfiles.keys():
		results = getfiles(search, test_data_dir)
		if len(results) == 1:
			demultfiles[search] = results[0]
		elif len(results) == 0:
			print("Error: Could not find file in", test_data_dir, "with", search, "in name")
			sys.exit()
		else:
			print("Found", len(results), "files with", search, "in name:")
			print("  - Using", results[0])
			demultfiles[search] = results[0]

	print("--------------------------------------------------------------")
	print(" Running Demultiplexor: Please wait for script to complete...")
	print("--------------------------------------------------------------")
	command = ["python","Demultiplexor.py","-r1",demultfiles["r1"],
				"-r2",demultfiles["r2"],
				"-i1",demultfiles["i1"],
				"-ix",demultfiles["ind"]]
	execute(command)
	return 1

def decombine():
	print ("-------------------------------------------------------------")
	print (" Running Decombinator: Please wait for script to complete...")
	print ("-------------------------------------------------------------")
	#find test fq files
	fqfiles = getfiles(["test","fq"], "./")
	print("Running for identified fastq files:")
	for f in fqfiles:
		print(os.path.basename(f))

	print("")
	for f in fqfiles:
	
		print("   Decombining", os.path.basename(f), "...")
		command = ["python","Decombinator.py",
				"-fq", os.path.basename(f)]
		execute(command)
		print("")
	return 1

def collapse():
	print("--------------------------------------------------------------")
	print(" Running Collapsinator: Please wait for script to complete...")
	print("--------------------------------------------------------------")
	#find test fq files
	fqfiles = getfiles(["dcr","n12"], "./")
	print("Running for identified dcr files:")
	for f in fqfiles:
		print(os.path.basename(f))

	print("")
	for f in fqfiles:
	
		print("   Collapsing", os.path.basename(f), "...")
		command = ["python","Collapsinator.py",
				"-in", os.path.basename(f), "-ol", "i8"]
		execute(command)
		print("")
	return 1

def translate():
	print("---------------------------------------------------------------")
	print(" Running CDR3Translator: Please wait for script to complete...")
	print("---------------------------------------------------------------")
	#find test fq files
	fqfiles = getfiles(["freq"], "./")
	print("Running for identified freq files:")
	for f in fqfiles:
		print(os.path.basename(f))

	print("")
	for f in fqfiles:
	
		print("   Translating", os.path.basename(f), "...")
		command = ["python","CDR3translator.py",
				"-in", os.path.basename(f)]
		execute(command)
		print("")
	return 1

if __name__ == '__main__':

	args = args()

	# test file directory
	test_data_dir = os.path.abspath("../Decombinator-Test-Data")
	if not os.path.isdir(test_data_dir):
		print("Error: could not find directory", test_data_dir)
		sys.exit()

	# check checkpoint is valid
	if args.checkpoint1.lower() not in ["demultiplexor", "decombinator", "collapsinator", "cdr3translator"]:
		print("Error: checkpoint1 argument must be one of demultiplexor/decombinator/collapsinator/cdr3translator")
		sys.exit()
	if args.checkpoint2.lower() not in ["demultiplexor", "decombinator", "collapsinator", "cdr3translator"]:
		print("Error: checkpoint2 argument must be one of demultiplexor/decombinator/collapsinator/cdr3translator")
		sys.exit()

	# run pipeline
	print("Pipeline will run scripts starting from", args.checkpoint1, "up to and including", args.checkpoint2)
	print("")

	# Demultiplexor
	if not args.checkpoint1.lower() in ["decombinator","collapsinator","cdr3translator"]:
		demultiplex(test_data_dir)

	if args.checkpoint2.lower() == "demultiplexor":
		sys.exit()

	# Decombinator
	if not args.checkpoint1.lower() in ["collapsinator","cdr3translator"]:
		decombine()
	if args.checkpoint2.lower() == "decombinator":
		sys.exit()

	# Collapsinator
	if not args.checkpoint1.lower() in ["cdr3translator"]:
		collapse()
	if args.checkpoint2.lower() == "collapsinator":
		sys.exit()

	# CDR3Translator
	translate()
	sys.exit()

		