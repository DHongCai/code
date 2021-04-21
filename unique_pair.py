#!/usr/bin/env python2.7
#--coding:utf-8--



##import package

import sys;
import os;
import re;

def switch_first(pgRNA,nc):
	if nc=="True":
		second_sg_loci=13
	else:
		second_sg_loci=15	
	sg1=pgRNA[0:second_sg_loci-1]
	sg2=pgRNA[(second_sg_loci-1):(second_sg_loci*2-2)]
	pgRNA[0:(second_sg_loci-1)]=sg2
	pgRNA[(second_sg_loci-1):(second_sg_loci*2-2)]=sg1
	return pgRNA


def select_bracode_pgRNA(pgRNAs,block_barcode,nc=False):
	barcode={}
	if nc=="True":
		second_sg_loci=13
	else:
		second_sg_loci=15

	pg_pool1=[]
	pg_pool2=[]
	pg_library={}
	##get all pgRNAs in to pool and library
	for i in range(0,len(pgRNAs)):
		pgRNA=pgRNAs[i].split("\n")[0].split("\t")
		pg_pool1+=[pgRNA]
		if pgRNA[1] not in pg_library:
			pg_library[pgRNA[1]]=0
		pg_library[pgRNA[1]]+=1
		if pgRNA[second_sg_loci] not in pg_library:
			pg_library[pgRNA[second_sg_loci]]=0
		pg_library[pgRNA[second_sg_loci]]+=1



	for pg in pg_pool1:
		sg1=pg[1]
		sg2=pg[second_sg_loci]
		if pg_library[sg1]==1 and sg1 not in barcode:
			if sg1 not in block_barcode:
				barcode[sg1]=pg
		elif pg_library[sg2]==1 and sg2 not in barcode:
			pg=switch_first(pg,nc)
			if sg2 not in block_barcode:
				barcode[sg2]=pg
		else:
			pg_pool2+=[pg]

	pg_library={}


	pg_pool1=pg_pool2
	pg_pool2=[]
	for pg in pg_pool1:
		sg1=pg[1]
		sg2=pg[second_sg_loci]

		if sg1 not in pg_library:
			pg_library[sg1]=0
		pg_library[sg1]+=1
	for pg in pg_pool1:
		sg1=pg[1]
		sg2=pg[second_sg_loci]
		if pg_library[sg1]==1 and sg1 not in barcode:
			if sg1 not in block_barcode:
				barcode[sg1]=pg

		elif sg2 not in pg_library and sg2 not in barcode:
			pg=switch_first(pg,nc)

			if sg2 not in block_barcode:
				barcode[sg2]=pg

		else:
			pg_pool2+=[pg]

	pg_library={}
	pg_pool1=pg_pool2
	pg_pool2=[]
	for pg in pg_pool1:
		sg1=pg[1]

		sg2=pg[second_sg_loci]		

		if sg1 not in barcode:
			if sg1 not in block_barcode:
				barcode[sg1]=pg

		elif sg2 not in barcode:
			pg=switch_first(pg,nc)
			if sg2 not in block_barcode:
				barcode[sg2]=pg

		else:
			pg_pool2+=[pg]



	# print ("Initial pool: "+str(len(pgRNAs))+", Seleted pool:"+str(len(barcode)))
	out_pgRNAs=[]
	for sg,pg in barcode.items():
		out_pgRNAs.append(pg)
	out_pgRNAs_sorted=sorted(out_pgRNAs,key = lambda x:float(x[-4]),reverse = True)
	#rerank the second different sgRNA to the top
	sg2_unique=[]
	out_pgRNAs_sorted_sg2sorted=[]
	for i in range(0,len(out_pgRNAs_sorted)):
		if out_pgRNAs_sorted[i][second_sg_loci] not in sg2_unique:
			sg2_unique.append(out_pgRNAs_sorted[i][second_sg_loci])
			if len(sg2_unique)>1:
				out_pgRNAs_sorted_sg2sorted.insert(len(sg2_unique),out_pgRNAs_sorted[i])
			else:
				out_pgRNAs_sorted_sg2sorted.append(out_pgRNAs_sorted[i])
		else:
			out_pgRNAs_sorted_sg2sorted.append(out_pgRNAs_sorted[i])




	

	return out_pgRNAs_sorted_sg2sorted

