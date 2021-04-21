"""
Process the data to count and normalize
"""
import argparse, tempfile, os, itertools, subprocess

def cut_adapter(file,outpath):
	outpath=(outpath+os.path.split(file)[0].split("/")[len(os.path.split(file)[0].split("/"))-1])
	mkdir(outpath)
	cmd="cutadapt --maximum-length=28 --minimum-length=17 --adapter={P1} --front={P2}  --output={P3} {P4}  ".format(
		P1="GTTTAAGAGCTATGCTGG",
		P2="TGGAAAGGACGAAACACCG",
		P3=outpath+"/"+os.path.split(file)[1].split(".")[0]+".fastq",
		P4=file
		)
	print cmd
	os.system(cmd)


def count(file,library,outpath):
	"""
	Read the fastq file and count the reads
	"""
	sgRNA2geneid_dict=sgRNA2geneid(library)



def sgRNA2geneid(library):
	"""
	Read sgRNA library 
	sgRNAid fa gene_ID
	"""
	genedict={}
	for line in open(library):
		ele=line.split("\n")[0].split("\t")
		genedict[ele[0]]=(ele[1],ele[2])
	return genedict


def GINI_index(x):
	xs=sorted(x)
	n=len(xs)
	gssum=sum([(i+1.0)*xs[i] for i in range(n)])
	ysum=sum(xs)
	if ysum==0.0:
		ysum=1.0
	gs=1.0-2.0*(n-gssum/ysum)/(n-1)
	return gs


def __main__():
	parser = argparse.ArgumentParser(
		description='{Process the fastq file from the CRISPR screening}')
	parser.add_argument('-o','--output',type=str,dest='outpath',help='output dir',default=os.getcwd())
	parser.add_argument('-i','--input_fastq',type=str,dest='input_fastq',help='The input fastq file')
	parser.add_argument('-l','--library',type=str,dest='library',help='The library file sgRNA_ID	sgRNA_fa	gene_ID')
	args=parser.parse_args()


	# cut_adapter()

	count(args.input_fastq,args.library,args.outpath)


"""
Commmand example
python 3_processData.py -i /picb/molsysbio/usr/caidh/projects/038_CD/WWSdata/3_removeAdapter/huh75_2/SRR4205827.fastq -l /picb/molsysbio/usr/caidh/projects/038_CD/WWSdata/4_summary/library.txt

"""

if __name__=="__main__":
	__main__()
