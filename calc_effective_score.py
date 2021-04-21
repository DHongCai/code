"""
The function to calculate the efficiency of sgRNA
"""

from sklearn.svm import SVC
import glob, copy, sys, math, operator, random, re, collections, tempfile, subprocess, logging, os, types
import pickle
from collections import defaultdict, Counter
from os.path import basename, join, splitext, isfile, dirname
import time, gzip, platform


# def calc_SSC_score(seqs):
# 	"""
# 	calc the effective score for each sgRNA
# 	"""
# 	continue

###Doench score	
def calc_doench_score(seq):
	guide_seq = seq[4:24]
	score = 0.59763615
	gcHigh = -0.1665878
	gcLow = -0.2026259
	gc_Num=guide_seq.count("G")+guide_seq.count("C")
	if gc_Num <=10:
		gc_Weight=gcLow
	if gc_Num>10:
		gc_Weight=gcHigh
	score+=abs(10-gc_Num)*gc_Weight

	for pos,modelSeq,weigth in Doench_weights:
		subSeq=seq[pos:pos+len(modelSeq)]
		if subSeq == modelSeq:
			score += weigth
	return 1.0/(1.0+math.exp(-score))

def calc_CrisprScan_score(seq):
	intercept=0.183930943629
	score=intercept
	for modelSeq,pos,weight in CRISPRscan_weights:
		subSeq=seq[pos-1:pos+len(modelSeq)-1]
		if subSeq==modelSeq:
			score+=weight
	return int(100*score)

def calc_Chari_score(seqs,outpath):
	Model_train="/picb/molsysbio/usr/caidh/projects/038_CD/scr/sgRNA.Scorer.1.0/293T.HiSeq.ST1.Nuclease.100.V2.SVM.Model.txt"
	seqs_index=Chari_seq_index(seqs)
	seqs_index_out=open(outpath+"seqs.index","w")
	seqs_index_out.write(seqs_index+"\n")
	seqs_index_out.close()
	cmd = "{P1} -v 0 {P2} {P3} {P4}".format(
		P1="/picb/molsysbio/usr/caidh/projects/038_CD/scr/Linux/svm_classify",
		P2=outpath+"seqs.index",
		P3=Model_train,
		P4=outpath+"seq.out")
	os.system(cmd)
	Chari_score=[]
	for line in open(outpath+"seq.out","r"):
		ele=line.split("\n")[0]
		Chari_score.append(ele)

	return Chari_score


def Chari_seq_index(seqs):
	vecs=[]
	for seq in seqs:
		vec=[]
		for i in range(0,len(seq)):
			for nu_index,char in enumerate("GTCA"):
				val=int(seq[i] ==char)
				if val != 0 :
					vec.append(("%d%d" %(i+1,nu_index+1),val))
		vecs.append(vec)
	lines=[]
	for vec in vecs:
		vec=["%s:%s" %(x,y) for x,y in vec]
		lines.append("0 "+" ".join(vec))
	return "\n".join(lines)







###datas
Doench_weights = [
# pasted/typed table from PDF and converted to zero-based positions
(1,'G',-0.2753771),(2,'A',-0.3238875),(2,'C',0.17212887),(3,'C',-0.1006662),
(4,'C',-0.2018029),(4,'G',0.24595663),(5,'A',0.03644004),(5,'C',0.09837684),
(6,'C',-0.7411813),(6,'G',-0.3932644),(11,'A',-0.466099),(14,'A',0.08537695),
(14,'C',-0.013814),(15,'A',0.27262051),(15,'C',-0.1190226),(15,'T',-0.2859442),
(16,'A',0.09745459),(16,'G',-0.1755462),(17,'C',-0.3457955),(17,'G',-0.6780964),
(18,'A',0.22508903),(18,'C',-0.5077941),(19,'G',-0.4173736),(19,'T',-0.054307),
(20,'G',0.37989937),(20,'T',-0.0907126),(21,'C',0.05782332),(21,'T',-0.5305673),
(22,'T',-0.8770074),(23,'C',-0.8762358),(23,'G',0.27891626),(23,'T',-0.4031022),
(24,'A',-0.0773007),(24,'C',0.28793562),(24,'T',-0.2216372),(27,'G',-0.6890167),
(27,'T',0.11787758),(28,'C',-0.1604453),(29,'G',0.38634258),(1,'GT',-0.6257787),
(4,'GC',0.30004332),(5,'AA',-0.8348362),(5,'TA',0.76062777),(6,'GG',-0.4908167),
(11,'GG',-1.5169074),(11,'TA',0.7092612),(11,'TC',0.49629861),(11,'TT',-0.5868739),
(12,'GG',-0.3345637),(13,'GA',0.76384993),(13,'GC',-0.5370252),(16,'TG',-0.7981461),
(18,'GG',-0.6668087),(18,'TC',0.35318325),(19,'CC',0.74807209),(19,'TG',-0.3672668),
(20,'AC',0.56820913),(20,'CG',0.32907207),(20,'GA',-0.8364568),(20,'GG',-0.7822076),
(21,'TC',-1.029693),(22,'CG',0.85619782),(22,'CT',-0.4632077),(23,'AA',-0.5794924),
(23,'AG',0.64907554),(24,'AG',-0.0773007),(24,'CG',0.28793562),(24,'TG',-0.2216372),
(26,'GT',0.11787758),(28,'GG',-0.69774)]


CRISPRscan_weights = [
# converted excel table of logistic regression weights with 1-based positions
('AA',18,-0.097377097),
('TT',18,-0.094424075),('TT',13,-0.08618771),('CT',26,-0.084264893),('GC',25,-0.073453609),
('T',21,-0.068730497),('TG',23,-0.066388075),('AG',23,-0.054338456),('G',30,-0.046315914),
('A',4,-0.042153521),('AG',34,-0.041935908),('GA',34,-0.037797707),('A',18,-0.033820432),
('C',25,-0.031648353),('C',31,-0.030715556),('G',1,-0.029693709),('C',16,-0.021638609),
('A',14,-0.018487229),('A',11,-0.018287292),('T',34,-0.017647692),('AA',10,-0.016905415),
('A',19,-0.015576499),('G',34,-0.014167123),('C',30,-0.013182733),('GA',31,-0.01227989),
('T',24,-0.011996172),('A',15,-0.010595296),('G',4,-0.005448869),('GG',9,-0.00157799),
('T',23,-0.001422243),('C',15,-0.000477727),('C',26,-0.000368973),('T',27,-0.000280845),
('A',31,0.00158975),('GT',18,0.002391744),('C',9,0.002449224),('GA',20,0.009740799),
('A',25,0.010506405),('A',12,0.011633235),('A',32,0.012435231),('T',22,0.013224035),
('C',20,0.015089514),('G',17,0.01549378),('G',18,0.016457816),('T',30,0.017263162),
('A',13,0.017628924),('G',19,0.017916844),('A',27,0.019126815),('G',11,0.020929039),
('TG',3,0.022949996),('GC',3,0.024681785),('G',14,0.025116714),('GG',10,0.026802158),
('G',12,0.027591138),('G',32,0.03071249),('A',22,0.031930909),('G',20,0.033957008),
('C',21,0.034262921),('TT',17,0.03492881),('T',13,0.035445171),('G',26,0.036146649),
('A',24,0.037466478),('C',22,0.03763162),('G',16,0.037970942),('GG',12,0.041883009),
('TG',18,0.045908991),('TG',31,0.048136812),('A',35,0.048596259),('G',15,0.051129717),
('C',24,0.052972314),('TG',15,0.053372822),('GT',11,0.053678436),('GC',9,0.054171402),
('CA',30,0.057759851),('GT',24,0.060952114),('G',13,0.061360905),('CA',24,0.06221937),
('AG',10,0.063717093),('G',10,0.067739182),('C',13,0.069495944),('GT',31,0.07342535),
('GG',13,0.074355848),('C',27,0.079933922),('G',27,0.085151052),('CC',21,0.088919601),
('CC',23,0.095072286),('G',22,0.10114438),('G',24,0.105488325),('GT',23,0.106718563),
('GG',25,0.111559441),('G',9,0.114600681)]