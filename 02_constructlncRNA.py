 #!/usr/bin/env python2.7
#--coding:utf-8 --


##import package
import glob,os,time,commands,shutil,sys
from joblib import Parallel,delayed
import re
from collections import defaultdict
import fileinput
import numpy as np
import HTSeq,pylab,pysam

##software
CUFFLINKS="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cufflinks"
GFFREAD="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/gffread"
CUFFCOMPARE="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cuffcompare"
CPAT="/picb/molsysbio/usr/caidh/packages/CPAT-1.2.2/bin/cpat.py"
CUFFMERGE="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cuffmerge"
CUFFQUANT="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cuffquant"
CUFFDIFF="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cuffdiff"
INTERSECT="/picb/molsysbio/apps/bin/bedtools/intersectBed"
RNAcode="/picb/molsysbio/usr/caidh/packages/RNAcode_install/bin/RNAcode"
PfamScan="/picb/molsysbio/usr/caidh/packages/Pfam/PfamScan/pfam_scan.pl"
TRANSEQ="/picb/molsysbio/usr/caidh/packages/EMBOSS-6.6.0/emboss/transeq"

##data
GTF="/picb/molsysbio2/caidh/035_Rmeta/8_conservation/data/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.97_chr.gtf"
GENOMEFA="/picb/molsysbio2/caidh/035_Rmeta/8_conservation/data/caenorhabditis_elegans/ce11.fa"
HEXA="/picb/molsysbio2/caidh/035_Rmeta/8_conservation/data/caenorhabditis_elegans/ce11_Hexamer.tsv"
Logit="/picb/molsysbio2/caidh/035_Rmeta/8_conservation/data/caenorhabditis_elegans/ce11.logit.RData"

##function
def mkdir(construct_outpath):	
	if not os.path.exists(construct_outpath):
		os.makedirs(construct_outpath)
def transcriptconstruct(file,outpath):
	sample=os.path.split(file)[1].split(".")[0]
	cmd="{P1}  -g {P2} --multi-read-correct  --seed 123 -p 12 -q -o {P3} {P4}".format(
		P1=CUFFLINKS,
		P2=GTF,
		P3=outpath+sample,
		P4=file
		)
	print cmd
	os.system(cmd)



def prepareGTFbed():
	fileout=open(GTF.replace("gtf","bed"),"w")
	ID=[]
	for line in open(GTF,"r"):
		if not re.search("^#",line):
			ele=line.split("\n")[0].split("\t")
			geneid=ele[8].split("gene_id \"")[1].split("\"")[0]
			if geneid not in ID:
				ID.append(geneid)
				genetype=ele[8].split("type \"")[1].split("\"")[0]
				lineout=[ele[0],ele[3],ele[4],geneid,genetype,ele[6]]
				fileout.write("\t".join(lineout)+"\n")
def getBed(bed):
	Ensemble=HTSeq.GenomicArrayOfSets("auto",stranded=False)
	#iv=HTSeq.GenomicInterval(ele[1],1,int(ele[3])+1000)
	for line in open(bed,"r"):
		ele=line.split("\n")[0].split("\t")
		start=int(ele[1])-500
		end=int(ele[2])+500
		# print line
		if start<0:
			start=1
		iv=HTSeq.GenomicInterval(ele[0],start,end)
		Ensemble[iv]+=ele[3]
	return Ensemble

def removeensemble(file,outpath):
	Ensemble=getBed(GTF.replace("gtf","bed"))	
	all_lncRNA=[]
	Non_Novel_lncRNA=[]
	gene_GTF=open(file,"r")
	for line in gene_GTF:
		ele=line.split("\n")[0].split("\t")
		# print ele
		
		if re.search("^CUFF",line) and float(ele[9])>0.1: ###FPKM more than 1 	and legnth >200
			# print ele[6]
			lncLen=int(ele[6].split(":")[1].split("-")[1])-int(ele[6].split(":")[1].split("-")[0])	
			if lncLen>200:
				all_lncRNA.append(ele[0])
				chrom=ele[6].split(":")[0]
				if not re.search("chr",chrom): ###
					chrom="chr"+chrom
				start=int(ele[6].split(":")[1].split("-")[0])
				end=int(ele[6].split(":")[1].split("-")[1])
				iv = HTSeq.GenomicInterval(chrom,start,end)
				# print iv
				overlap=set()
				# print line
				for key,value in Ensemble[iv].steps():
					overlap.update(value)

				if len(overlap)>0:
					# print overlap
					# print ele[0]
					Non_Novel_lncRNA.append(ele[0])

	Novel_lncRNA=list(set(all_lncRNA).difference(set(Non_Novel_lncRNA)))

	gene_GTF.close()
	print "NovellncRNA "+str(len(Novel_lncRNA))


	transcriptID2exonNum={}
	trans_GTF=open(os.path.split(file)[0]+"/transcripts.gtf","r")
	for line in trans_GTF:
		ele=line.split("\n")[0].split("\t")
		if ele[2] =="exon":
			geneid=ele[8].split("gene_id \"")[1].split("\"")[0]
			transcriptID=ele[8].split("transcript_id \"")[1].split("\"")[0]

			# FPKM=float(ele[8].split("FPKM \"")[1].split("\"")[0])
			if geneid in Novel_lncRNA:
				# print geneid
				# print transcriptID
				if transcriptID not in transcriptID2exonNum:
					transcriptID2exonNum[transcriptID]=1
				else:
					transcriptID2exonNum[transcriptID]=transcriptID2exonNum[transcriptID]+1
	Geneid_singleexon=[]
	Geneid_multiexon=[]
	trans_GTF.close()
	for key,value in transcriptID2exonNum.items():
		# print key
		if value==1:
			key0=key.split(".")
			if len(key0)>1:
				geneid=key0[0]+"."+key0[1]
				Geneid_singleexon.append(geneid)
		else:
			key0=key.split(".")
			if len(key0)>1:	
				geneid=key0[0]+"."+key0[1]
				Geneid_multiexon.append(geneid)
	Geneid_multiexon=list(set(Geneid_multiexon).difference(set(Geneid_singleexon)))			

	print "singleexon "+str(len(Geneid_singleexon))
	print "multiexon "+str(len(Geneid_multiexon))
	fileout=open(os.path.split(file)[0]+"/transcripts_filterOverlap.gtf","w")
	trans_GTF=open(os.path.split(file)[0]+"/transcripts.gtf","r")
	for line in trans_GTF:
		ele=line.split("\n")[0].split("\t")
		if not re.search("chr",ele[0]):
			line="chr"+line
		geneid=ele[8].split("gene_id \"")[1].split("\"")[0]
		FPKM=float(ele[8].split("FPKM \"")[1].split("\"")[0])
		if geneid in Geneid_singleexon and FPKM>5:
			fileout.write(line.split("\n")[0]+" exon \"single\";\n")
		if geneid in Geneid_multiexon:
			fileout.write(line.split("\n")[0]+" exon \"multi\";\n")
	trans_GTF.close()


def callCPAT(file,outpath):
	fileout=file.split(".")[0]+".fa"
	cmd="{P1} -w {P2} -g {P3} {P4} ".format(
		P1=GFFREAD,
		P2=fileout,
		P3=GENOMEFA,
		P4=file
		)
	print cmd
	os.system(cmd)
	cmd="{P1} -g {P2} -d {P3} -x {P4} -o {P5}".format(
		P1=CPAT,
		P2=file.split(".")[0]+".fa",
		P3=Logit,
		P4=HEXA,
		P5=file.split(".")[0]+".cpat"
		)
	print cmd
	os.system(cmd)
	filein =open(file.split(".")[0]+".cpat","r")
	next(filein)
	Hexacutoff=0
	codingcutoff=0.44
	codingpotentiallnc={}
	for line in filein:

		ele=line.split("\n")[0].split("\t")
		if  float(ele[5]) < 0.44 and float(ele[4]) < 0:
			codingpotentiallnc[ele[0]]=1
			#print ele[0] 
	fileout=open(file.split(".")[0]+"_CPATfilter.gtf","w")
	print file.split(".")[0]+"_CPATfilter.gtf"
	for line in open(file.split(".")[0]+".gtf","r"):
		ele=line.split("\n")[0].split("\t")
		attr=ele[8].split(" ")
		transcriptid=attr[3].replace("\"","").replace(";","")
		if transcriptid in codingpotentiallnc:
			fileout.write(line)


def callPfamscan(file,outpath):
	cmd="{P1} -trim {P2} {P3}".format(P1=TRANSEQ,P2=file.split(".")[0]+".fa",P3=file.split(".")[0]+"_protein.fa")
	print cmd
	os.system(cmd)
	# cmd="sed 's/*//g' {P1} >{P2} |rm {P1}".format(P1=file.split(".")[0]+"_Midprotein.fa",P2=file.split(".")[0]+"_protein.fa")
	# print cmd
	# os.system(cmd)
	cmd="{P1}  -cpu 20 -fasta {P2}  -dir  /picb/molsysbio/usr/caidh/packages/Pfam/PfamScan/data  -out {P3}".format(
		P1=PfamScan,
		P2=file.split(".")[0]+"_protein.fa",
		P3=file.split(".")[0]+"_Pfamfilter.txt")
	print cmd
	os.system(cmd)
	codinglnc=[]
	for line in open(file.split(".")[0]+"_Pfamfilter.txt","r"):
		if re.search("^CUFF",line):
			ele=line.split("\n")[0].split("\t")
			lnc=ele[0].split(".")[0]+"."+ele[0].split(".")[1]
			if lnc not in codinglnc:
				codinglnc.append(lnc)

	fileout=open(file.split(".")[0]+"_CPATfilter_Pfamfilter.gtf","w")
	for line in open(file.split(".")[0]+"_CPATfilter.gtf","r"):
		ele=line.split("\n")[0].split("\t")
		attr=ele[8].split(" ")
		geneid=attr[1].replace("\"","").replace(";","")
		if geneid not in codinglnc:
			fileout.write(line)


"""
gffread -w ./output/gffcompare/lncRNA_transcript.fa -g ./genome/GCF_000003025.6_Sscrofa11.1_genomic.fna ./output/gffcompare/assembly_lncNRA.gtf
transeq ./output/gffcompare/lncRNA_transcript.fa ./output/gffcompare/lncRNA_protein_mid.fa
sed 's/*//g' ./output/gffcompare/lncRNA_protein_mid.fa > ./output/gffcompare/lncRNA_protein.fa
rm ./output/gffcompare/lncRNA_protein_mid.fa
./pfam_scan.pl -cpu 20 -fasta ./output/gffcompare/lncRNA_protein.fa -dir ~/disk_4T/Pfam_database -out ./output/PfamScan/Pfam_result.txt
less Pfam_result.txt|grep -v '#'|awk '{print $1}' > Pfam_id.txt
comm -13 <(sort Pfam_id.txt) <(sort lncRNA_id.txt) > Pfam_id1.txt
mv Pfam_id1.txt Pfam_id.txt
"""

def makecuffmerge(files,outpath):
	mkdir(outpath+"4_novolnc/")
	fileout=open(outpath+"4_novolnc/lnc_merge.list","w")
	for file in files:
		fileout.write(file+"\n")

	cmd="{P1} --keep-tmp  -o {P3} -p 5 {P2} ".format(
		P1=CUFFMERGE,
		P2=outpath+"4_novolnc/lnc_merge.list",
		P3=outpath+"4_novolnc/novo"
		)
	print cmd
	os.system(cmd)


def combine_GTF(novoGTF,index):
	G2T={}
	G2G={}
	fileout=open(novoGTF.replace(".gtf","_new.gtf"),"w")
	for line in open(novoGTF,"r"):
		line=line.split("\n")[0]	
		line=line.replace("CUFF.",index+"_")
		line=line+" gene_biotype \"lncRNAnovo\";"
		ele=line.split("\n")[0].split("\t")
		ele[8]=ele[8].replace(".","_")
		line="\t".join(ele)
		geneid=ele[8].split("gene_id \"")[1].split("\"")[0]

		if not re.search("_",ele[0]):
			if geneid not in G2T:
				G2T[geneid]=[]
				G2G[geneid]=ele
				G2G[geneid][2]="gene"
			G2T[geneid].append(line)
			if int(ele[3])<	int(G2G[geneid][3]):
				G2G[geneid][3]=ele[3]
			if int(ele[4])>	int(G2G[geneid][4]):
				G2G[geneid][4]=ele[4]
	G2G_key=sorted(G2G)
	for key in G2G_key:
		fileout.write("\t".join(G2G[key])+"\n")
		for line in G2T[key]:
			fileout.write(line+"\n")

	cmd="cat {P1} {P2} >{P3}".format(
		P1=GTF,
		P2=novoGTF.replace(".gtf","_new.gtf"),
		P3=GTF.replace(".gtf","_addNovo.gtf"))
	print cmd
	# os.system(cmd)

def HTSEQcount(file,outpath):
	cmd_sense_exon="htseq-count --stranded=no -f bam --order=pos --type=exon --idattr=gene_id {P2} {P3} > {P4}".format(
		P2=file,
		P3=GTF.replace(".gtf","_addNovo.gtf"),
		P4=outpath+os.path.split(file)[1].split(".")[0].replace("Aligned","")+".sense_exon")
	# print cmd_sense_exon
	# os.system(cmd_sense_exon)
	cmd_sense_gene="htseq-count --stranded=no -f bam --order=pos --type=gene --idattr=gene_id {P2} {P3} > {P4}".format(
		P2=file,
		P3=GTF.replace(".gtf","_addNovo.gtf"),
		P4=outpath+os.path.split(file)[1].split(".")[0].replace("Aligned","")+".sense_gene")
	# print cmd_sense_gene
	# os.system(cmd_sense_gene)
	cmd="stringtie {P2} -G {P3} -p 20 -A {P4} -o {P5}".format(
		P2=file,
		P3=GTF.replace(".gtf","_addNovo_nochr.gtf"),
		P4=outpath+os.path.split(file)[1].split(".")[0].replace("Aligned","")+".count",
		P5=outpath+os.path.split(file)[1].split(".")[0].replace("Aligned","")+".gtf")
	print cmd
	os.system(cmd)
	
import pandas as pd	
def makeintomatrix(files,outpath):
	datas={}

	for file in files:
		data={}
		name=os.path.split(file)[1].split(".")[0]
		print file
		i=0
		for line in open(file,"r"):
			ele=line.split("\n")[0].split("\t")
			gene=ele[0]+"|"+ele[2]+":"+ele[4]+"-"+ele[5]
			if i != 0 :
				if not re.search("STRG",line):
					data[gene]=ele[8]
			i=1
		datas[name]=data
	matrix=pd.DataFrame(datas)
	matrix.fillna(0)
	matrix.to_csv(outpath+"Allgene.expression",sep="\t",index_label="gene")


def main():
	outpath="/picb/molsysbio2/caidh/032_ASL/GSE101964_celegans_glp1_riboM/"
	# mkdir(outpath+"3_construct_lncRNA_2/")
	# bams=glob.glob(outpath+"2_mapping/*/*bam")
	# Parallel(n_jobs=len(bams))(delayed(transcriptconstruct)(bam,outpath+"3_construct_lncRNA_2/") for bam in bams)

	# prepareGTFbed()

	# files=glob.glob(outpath+"3_construct_lncRNA_2/*/genes.fpkm_tracking")
	# Parallel(n_jobs=len(files))(delayed(removeensemble)(file,outpath)for file in files)


	# ## Coding pontential
	# files=glob.glob(outpath+"3_construct_lncRNA_2/*/*filterOverlap.gtf")
	# Parallel(n_jobs=10)(delayed(callCPAT)(file,outpath)for file in files)

	#Pfamscan
	# files=glob.glob(outpath+"3_construct_lncRNA_2/*/*filterOverlap.gtf")
	# Parallel(n_jobs=20)(delayed(callPfamscan)(file,outpath)for file in files)

	#CUFFMerge
	# files=glob.glob(outpath+"3_construct_lncRNA_2/*/*_CPATfilter_Pfamfilter.gtf")
	# makecuffmerge(files,outpath)


	# combine_GTF(outpath+"4_novolnc/novo/transcripts.gtf","celegans")
	# mkdir(outpath+"5_HTseq/")
	bams=glob.glob(outpath+"2_mapping/*/*bam")
	# bams.sort()

	Parallel(n_jobs=10)(delayed(HTSEQcount)(bam,outpath+"5_HTseq/") for bam in bams)	



	##make expression matrix
	mkdir(outpath+"6_summary/")
	GTFs=glob.glob(outpath+"5_HTseq/*count")
	makeintomatrix(GTFs,outpath+"6_summary/")



















if __name__=="__main__":
	main()

