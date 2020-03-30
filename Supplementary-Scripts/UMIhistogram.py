import sys
import os
import matplotlib.pyplot as plt
import argparse

def args():
	parser = argparse.ArgumentParser(description='Histogram plotting tool for average UMI cluster sizes, post collapsing')
	parser.add_argument('-in', '--infile', type=str, help='File containing histogram data', required=True)
	parser.add_argument('-o', '--outfile', type=str, help='Ouput png file name (default automatically names file)', required=False)
	parser.add_argument('-b', '--bins', type=int, help='Number of equal width bins', required=False, default=100)
	parser.add_argument('-c', '--color', type=str, help='Colour of histogram bars', required=False, default=None)	
	parser.add_argument('-d', '--dpi', type=int, help='DPI of output figure', required=False, default=None)	
	return parser.parse_args()

def checkInfile(infile):
	if not os.path.exists(args.infile):
		print("Error: could not find", infile)
		print("Please check input filename and which directory you are in.")
		sys.exit()
	return 1

def getOutfileName(infile):
	# find where Decombinator repository is stored
	cwd = os.getcwd()
	outdir = cwd.split("Decombinator")[0] + os.sep + "Decombinator" + os.sep + "Logs"
    
	# Logs directory should exist, but if not...
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	outfile = os.path.splitext( os.path.basename(infile) )[0] + ".png"

	return outdir + os.sep + outfile

def getData(infile):
	av_cluster_sizes = []
	counts = []
	with open(infile, "rt") as f:
		for line in f:
			clus_size, c = line.split(",")
			av_cluster_sizes.append(int(clus_size))
			counts.append(int(c))
	return (av_cluster_sizes, counts)

def histogram(data, outfilename, bins = 100, color = None, dpi = None):
	plt.hist(data[0], weights = data[1], bins = bins, color = color)
	plt.xlabel("Average UMI cluster size")
	plt.ylabel("Frequency")
	plt.savefig(outfilename, dpi = dpi)
	print("Histogram saved to", outfilename)

if __name__ == '__main__':

	args = args()
	checkInfile(args.infile)
	
	if args.outfile:
		outfilename = args.outfile
	else:
		outfilename = getOutfileName(args.infile)
	
	data = getData(args.infile)

	histogram = histogram(data, outfilename, args.bins, args.color, args.dpi)
	