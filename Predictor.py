"""
Predict the effeciency of pgRNA

"""


from __future__ import division
import sys
import argparse
import os

from collections import defaultdict
from sklearn.svm import SVC
import pickle
import warnings

good_pg="/picb/molsysbio/usr/caidh/projects/038_CD/WWSdata/5_pgRNA_feature/3_effective_sgRNA/good_pgRNA.list"
bad_pg="/picb/molsysbio/usr/caidh/projects/038_CD/WWSdata/5_pgRNA_feature/3_effective_sgRNA/bad_pgRNA.list"
outpath="/picb/molsysbio/usr/caidh/projects/038_CD/WWSdata/5_pgRNA_feature/3_effective_sgRNA/"
#binary encoding
encoding= defaultdict(str)
encoding['A']='0001'
encoding['C']='0010'
encoding['T']='0100'
encoding['G']='1000'


def learning(good,bad):
	xList=[]
	yList=[]
	for line in open(good,"r"):
		ele=line.split("\n")[0].split("\t")
		if int(ele[15])<int(ele[20]):
			sequence1=ele[17]
			sequence2=ele[22]
		else:
			sequence1=ele[22]
			sequence2=ele[17]		
		sequence=sequence1+sequence2
		xList.append(index_sequence(sequence))
		yList.append(1)
	for line in open(bad,"r"):
		ele=line.split("\n")[0].split("\t")
		if int(ele[15])<int(ele[20]):
			sequence1=ele[17]
			sequence2=ele[22]
		else:
			sequence1=ele[22]
			sequence2=ele[17]	
		sequence=sequence1+sequence2
		xList.append(index_sequence(sequence))
		yList.append(-1)
	
	#learning
	SVM_learn=SVC(kernel='linear')
	SVM_learn.fit(xList,yList)

	return(SVM_learn)

def index_sequence(seq):
	index=[]
	for x in range(0,40):
		for y in range(0,4):	
			index.append(int(encoding[seq[x]][y]))
	return index


def predict(sequence):
	warnings.filterwarnings("ignore")
	SVM_learn = pickle.load(open(outpath+"learning_model.pickle", 'rb'))
	seq_index=index_sequence(sequence)
	score= SVM_learn.decision_function(seq_index)
	return float(score)



def main():
	outpath="/picb/molsysbio/usr/caidh/projects/038_CD/WWSdata/5_pgRNA_feature/3_effective_sgRNA/"
	# SVM_learn=learning(good_pg,bad_pg)
	# pickle.dump(SVM_learn,open(outpath+"learning_model.pickle","wb"))

	# predict("GGTGAACGGCTGCGCGACAGGGTGAACGGCTGCGCGACAG")



if __name__=="__main__":
	main()
