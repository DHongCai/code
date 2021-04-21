#!/usr/bin/env python2.7
#--coding:utf-8--



##import package
import glob,os,time,commands,shutil,sys
from joblib import Parallel,delayed
import argparse, tempfile, os, itertools, subprocess
from Predictor import *
import HTSeq
def select_good_sgRNA(sgRNAs,noncoding):
	sgRNAs_out=[]
	coding=get_coding_region()
	i=0
	for line in open(sgRNAs,"r"):
		if i>0:
			ele=line.split("\n")[0].split("\t")
			if pick_nucletide(ele[1]) and float(ele[7])>10 and float(ele[8])>0.1:## efficiency and the off target value 

				if noncoding=="y":##remove coding region
					iv=HTSeq.GenomicInterval(ele[3],int(ele[6]),int(ele[6])+1,".")
					
					for iv1, gene in coding[iv].steps():
						if len(gene)==0:
							sgRNAs_out.append("_".join(ele))
				else:
					sgRNAs_out.append("_".join(ele))			
		i=i+1
	selected_num=len(sgRNAs_out)
	return sgRNAs_out,i,selected_num



def pick_nucletide(guide):
	"""
	return true if guide does not contain any of the following homopolymer target strech
	"""
	guide=guide.upper()
	# print guide
	st1 = ('AAAA')
	st2 = ('TTTT')
	st3 = ('GGGG')
	st4 = ('CCCC')
	st5 = ('CGTCTC')
	st6 = ('GAGACG')
	if not ((st1 in guide) or (st2 in guide) or (st3 in guide) or (st4 in guide) or (st5 in guide) or (st6 in guide)):
		return 'true'



def make_allpairs(guides):
	allpairs=[]
	gtf_index=get_transtript_location()
	number_beforeTSS=get_before_TSS(guides,gtf_index)
	No_cross_TSS_num=0

	for i in range(0,len(guides)):
		sgi=guides[i].split("_")
		for j in range(i+1,len(guides)):
			sgj=guides[j].split("_")
			
			gene_gtf=gtf_index[sgi[0].split("|")[1]]
			if gene_gtf[5]=="+":
				tss=int(gene_gtf[1])
			else:
				tss=int(gene_gtf[2])
			if number_beforeTSS[0]>5 and number_beforeTSS[1]>5:
				if (int(sgi[6]) > tss and int(sgj[6])<tss) or (int(sgi[6]) < tss and int(sgj[6])>tss):
					allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
			elif number_beforeTSS[1]==0:
				if No_cross_TSS_num<1000:
					if abs(int(sgi[6])-int(sgj[6]))>1000 and abs(int(sgi[6])-int(sgj[6]))<10000 and int(sgi[6])<int(sgj[6]):
						allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
						No_cross_TSS_num=No_cross_TSS_num+1

			elif number_beforeTSS[0]==0:
				if No_cross_TSS_num<1000:
					allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
					No_cross_TSS_num=No_cross_TSS_num+1
			else:
				if (int(sgi[6]) > tss and int(sgj[6])<tss) or (int(sgi[6]) < tss and int(sgj[6])>tss):
					allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
				if No_cross_TSS_num<1000:
					if gene_gtf[5]=="+" and (int(sgi[6]) > tss and int(sgj[6])>tss) and int(sgi[6])<int(sgj[6]):
						allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
					elif gene_gtf[5]=="-" and (int(sgi[6]) < tss and int(sgj[6])<tss) and int(sgi[6])<int(sgj[6]):
						allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
					elif abs(int(sgi[6])-int(sgj[6]))>1000 and abs(int(sgi[6])-int(sgj[6]))<10000 and int(sgi[6])<int(sgj[6]):
						allpairs.append("_".join(sgi)+"&&"+"_".join(sgj))
						No_cross_TSS_num=No_cross_TSS_num+1

					No_cross_TSS_num=No_cross_TSS_num+1
	return allpairs,number_beforeTSS

def get_before_TSS(guides,gtf_index):
	num=[0,0]
	for i in range(0,len(guides)):
		sgi=guides[i].split("_")
		gene_gtf=gtf_index[sgi[0].split("|")[1]]
		if gene_gtf[5]=="+":
			tss=int(gene_gtf[1])
			if int(sgi[6]) < tss:
				num[0]+=1
			else:
				num[1]+=1
		else:
			tss=int(gene_gtf[2])
			if int(sgi[6]) > tss:
				num[0]+=1
			else:
				num[1]+=1
	return num



def give_efficiency(pairs):
	all_pairs_effi=[]
	for i in range(0,len(pairs)):
		all_pairs=pairs[i].split("\n")[0].split("&&")
		all_pairs1=all_pairs[0].split("_")
		all_pairs2=all_pairs[1].split("_")
		distance=abs(int(all_pairs1[6])-int(all_pairs2[6]))
		if int(all_pairs1[6])<int(all_pairs2[6]):
			pg=all_pairs1[1]+all_pairs2[1]
			# score=predict(pg)
			score=0
		else:
			pg=all_pairs2[1]+all_pairs1[1]
			# score=predict(pg)		
			score=0

		all_pairs_effi.append(pairs[i].split("\n")[0]+"||"+str(distance)+"_"+str(score)+"\n")
	return all_pairs_effi

def give_feature(pairs):
	all_pairs_effi_feature=[]

	cds=get_cds_region()
	exon=get_exon_region()
	gtf_index=get_transtript_location()
	for i in range(0,len(pairs)):
		first_feature=''
		second_feature=''
		both_feature=''
		all_pairs=pairs[i].split("\n")[0].split("&&")
		all_pairs1=all_pairs[0].split("_")
		sgRNA1_iv=HTSeq.GenomicInterval(all_pairs1[3],int(all_pairs1[6]),int(all_pairs1[6])+1,".")
		all_pairs2=all_pairs[1].split("_")	
		sgRNA2_iv=HTSeq.GenomicInterval(all_pairs2[3],int(all_pairs2[6]),int(all_pairs2[6])+1,".")
		distance=abs(int(all_pairs1[6])-int(all_pairs2[6]))
		gene_gtf=gtf_index[all_pairs1[0].split("|")[1]]
		strand=gene_gtf[5]
		if all_pairs1[0].split("|")[1] in cds:
			cds_start=int(cds[all_pairs1[0].split("|")[1]])
		else:
			if strand=="+":
				cds_start=gene_gtf[2]
			else:
				cds_start=gene_gtf[1]


		if (strand=="+" and all_pairs1[2]=="t" and all_pairs2[2]=="t") or (strand=="-" and all_pairs1[2]=="b" and all_pairs2[2]=="b"):
			both_feature+='ORI;'
		else:
			both_feature+='NOTORI;'
		if distance<10000 and distance>500:
			both_feature+='GOODDIS;'
		elif distance<501:
			both_feature+='CLOSEDIS;'
		else:
			both_feature+='FARDIS;'

		if strand=="+":
			tss=int(gene_gtf[1])
			if int(all_pairs1[6])<tss:
				first_feature+='TSS;'
			else:
				first_feature+='NOTSS;'
			if int(all_pairs2[6])<tss:
				second_feature+='TSS;'
			else:
				second_feature+='NOTSS;'
		else:
			tss=int(gene_gtf[2])
			if int(all_pairs1[6])>tss:
				first_feature+='TSS;'
			else:
				first_feature+='NOTSS;'
			if int(all_pairs2[6])>tss:
				second_feature+='TSS;'
			else:
				second_feature+='NOTSS;'			
				
		if strand=="+":
			if int(all_pairs1[6])>cds_start:
				first_feature+='CDS;'
			else:
				first_feature+='NOTCDS;'
			if int(all_pairs2[6])>cds_start:
				second_feature+='CDS;'
			else:
				second_feature+='NOTCDS;'
		else:
			if int(all_pairs1[6])<cds_start:
				first_feature+='CDS;'
			else:
				first_feature+='NOTCDS;'
			if int(all_pairs2[6])<cds_start:
				second_feature+='CDS;'
			else:
				second_feature+='NOTCDS;'			


		
		for iv1, gene in exon[sgRNA1_iv].steps():
			if len(gene)>0:
				first_feature+='EXON;'
			
			else:
				first_feature+='NOEXON;'
			break
		
		for iv1, gene in exon[sgRNA2_iv].steps():
			if len(gene)>0:
				second_feature+='EXON;'

			else:
				second_feature+='NOEXON;'
			break
		all_pairs_effi_feature.append(pairs[i].split("\n")[0]+"_"+first_feature+"_"+second_feature+"_"+both_feature+"\n")
	return(all_pairs_effi_feature)




def get_transtript_location():
	gtf="/picb/molsysbio/usr/caidh/data/hg38/gencode.v24.annotation.bed"
	gtf_index={}
	for line in open(gtf,"r"):
		ele=line.split("\n")[0].split("\t")
		gtf_index[ele[3].split("|")[1]]=ele
	return gtf_index

def get_coding_region():
	bed="/picb/molsysbio/usr/caidh/data/hg38/gencode_coding.v24.annotation.bed"
	coding=HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for line in open(bed,"r"):
		ele=line.split("\n")[0].split("\t")
		if ele[5]=="+":
			iv=HTSeq.GenomicInterval(ele[0],int(ele[1])-5000,int(ele[2]),".")
			coding[iv]+=ele[3]
		else:
			iv=HTSeq.GenomicInterval(ele[0],int(ele[1]),int(ele[2])+5000,".")
			coding[iv]+=ele[3]
	return coding 



def get_cds_region():
	bed="/picb/molsysbio/usr/caidh/data/hg38/gencode_cds.v24.annotation.bed"
	cds={}
	for line in open(bed,"r"):
		ele=line.split("\n")[0].split("\t")

		if ele[5]=="+":
			if ele[3] not in cds:
				cds[ele[3]]=ele[1]
			else:
				if int(ele[1])<int(cds[ele[3]]):
					cds[ele[3]]=ele[1]
		if ele[5] =="-":
			if ele[3] not in cds:
				cds[ele[3]]=ele[2]
			else:
				if int(ele[2])>int(cds[ele[3]]):
					cds[ele[3]]=ele[2]
	return cds
def get_exon_region():
	bed="/picb/molsysbio/usr/caidh/data/hg38/gencode_keepExon.v24.annotation.bed"
	exon=HTSeq.GenomicArrayOfSets("auto",stranded=False)
	for line in open(bed,"r"):
		ele=line.split("\n")[0].split("\t")
		iv=HTSeq.GenomicInterval(ele[0],int(ele[1]),(int(ele[2])+1),".")
		exon[iv]+=ele[3]
	return exon 	
def writesgRNA(all_pairs_effi_feature,outpath):
	outfile=open(outpath+"/pgRNA.design","w")

	outfile.write("\t".join(["geneid","sgRNA","strand","chrom","start","end","location","off_target_index","SSC_effective_value","Doench_score","CRISPRscan_score","Chari_score",
		"geneid","sgRNA","strand","chrom","start","end","location","off_target_index","SSC_effective_value","Doench_score","CRISPRscan_score","Chari_score","Distance","pgRNA_efficiency","first_feature","second_feature","all_feature"])+"\n")
	if len(all_pairs_effi_feature)>1:
	# 	all_pairs_effi_sorted=sorted(all_pairs_effi_feature,key = lambda x:float(x.split("_")[-4]),reverse = True)
		for i in range(0,len(all_pairs_effi_feature)):
			lineout=all_pairs_effi_feature[i].replace("_","\t").replace("&&","\t").replace("||","\t")
			outfile.write(lineout)



def __main__():
	parser = argparse.ArgumentParser(
		description='Make sgRNA intop pgRNA')
	parser.add_argument('-o','--output',type=str,dest='outpath',help='output dir',default=os.getcwd())
	parser.add_argument('-g','--sgRNAs',type=str,dest='sgRNAs',help='the list of sgRNAs')
	parser.add_argument('-n','--noncoding',type=str,dest='noncoding',help='noncoding region: y or n')

	args=parser.parse_args()

	log=open(args.outpath+"/design_pgRNA.log","w")
	guides,previous_num,selected_num = select_good_sgRNA(args.sgRNAs,args.noncoding)
	log.write("all sgRNA: "+str(previous_num)+"\n")
	log.write("selected sgRNA: "+str(selected_num)+"\n")
	print "all sgRNA: "+str(previous_num)
	print "selected sgRNA: "+str(selected_num)	

	if selected_num>5:
		log.write("Begin to make pairs gRNA."+"\n")
		all_pairs,number_beforeTSS=make_allpairs(guides)
		log.write("sgRNA ahead TSS: "+str(number_beforeTSS[0])+"\n")
		log.write("sgRNA after TSS: "+str(number_beforeTSS[1])+"\n")
		log.write("Design pairs number: "+str(len(all_pairs))+"\n")
		print "sgRNA ahead TSS: "+str(number_beforeTSS[0])
		print "sgRNA after TSS: "+str(number_beforeTSS[1])
		print "Design pairs number: "+str(len(all_pairs))

		all_pairs_effi=give_efficiency(all_pairs)

		all_pairs_effi_feature=give_feature(all_pairs_effi)

		writesgRNA(all_pairs_effi_feature,args.outpath)
	else:
		log.write("no suitable sgRNA selected to make pairs."+"\n")

"""

python 2_pgRNA_design.py -g sgRNA.design -n y



"""

if __name__=="__main__":
	__main__()

