#!/usr/bin/env python2.7
#--coding:utf-8 --


##import package
import glob,os,time,commands,shutil,sys
from joblib import Parallel,delayed
import re

##software
TOPHAT2="/picb/molsysbio/usr/xuanhongwen/tool/tophat-2.0.13.Linux_x86_64/tophat2"
STAR="/picb/molsysbio/usr/caidh/packages/STAR-master/bin/Linux_x86_64_static/STAR"
cutadapter1="/picb/molsysbio/usr/caidh/packages/SG/sage-6.2/local/bin/cutadapt"
CUFFLINKS="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cufflinks"
CUFFDIFF="/picb/molsysbio/usr/chenghao/SoftWare/cufflinks-2.2.1.Linux_x86_64/cuffdiff"
BAM2WIG="bam2wig.py"
WIG2BIGWIG="/picb/molsysbio/usr/caidh/packages/wigToBigWig"
SAMTOOLS="/picb/molsysbio/usr/chenghao/SoftWare/samtools-1.1/samtools"
STAR_INDEX="/picb/molsysbio/usr/caidh/data/ce/STAR"
##file
inpath="/picb/molsysbio/usr/xuanhongwen/set-32/rnaseq/data/0.fastq/"
outpath="/picb/molsysbio/usr/caidh/projects/017_celegans/result/1_processdata/"
fastq_file="2_fastq/"
mapping="3_mapping/"
stat_bed_bigwig="4_stat_bed_bigwig/"
count="5_count/"
DEG="6_DEG"
BOWTIE_INDEX="/picb/molsysbio/usr/caidh/data/ce/bowtieindex/WS254"
GTF="/picb/molsysbio/usr/caidh/projects/017_celegans/data/c_elegans.PRJNA13758.WS254.canonical_geneset.gtf"
FA="/picb/molsysbio/usr/caidh/data/ce/c_elegans.PRJNA13758.WS254.genomic_softmasked.fa"
CE10SIZE="/picb/molsysbio/usr/caidh/data/chrom.size/ce10.chrom.sizes"
##function
def prepare_files(path):
	for p in [path+fastq_file,path+mapping,path+stat_bed_bigwig,path+count,path+DEG]:
		if not os.path.exists(p):
			os.mkdir(p)
def sra2fastq(file,outpath):
	mkdir(outpath+os.path.split(file)[1].split(".")[0])
	cmd="fastq-dump --gzip --split-3  {P1} -O {P2}".format(P1=file,P2=outpath+os.path.split(file)[1].split(".")[0])
	print cmd
	os.system(cmd)

def runFastqc(fq):
    dir = os.path.split(fq)[0]
    #print dir
    odir = '%s/Fastqc' % (dir)
    # if not os.path.exists(odir):
    #     #return
    # 	os.mkdir(odir)
    cmd = 'fastqc -o %s -f fastq %s' % (odir,fq)
    print cmd
    os.system(cmd)

def mkdir(construct_outpath):	
	if not os.path.exists(construct_outpath):
		os.makedirs(construct_outpath)
def STARmappingsingle(file,pathin):
	# print file
	"""
STAR --runMode alignReads --runThreadN 10 --genomeDir /picb/molsysbio/usr/chenghao/GeneAnnotation/STAR/mm10 --readFilesIn ../RawData/P15R451.fastq --sjdbGTFfile /picb/molsysbio/usr/chenghao/GeneAnnotation/mm10/mm10.refgene.gt
f --clip5pNbases 10 --outFileNamePrefix DICER_WT-0h --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outWigStrand Stranded --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --quantMode GeneCounts 

STAR   --runMode genomeGenerate   --runThreadN 5   --genomeDir /picb/molsysbio/usr/chenghao/GeneAnnotation/STAR/mm10   --genomeFastaFiles /picb/molsysbio/usr/chenghao/GeneAnnotation/mm10/mm10.fa   

--sjdbGTFfile {P4} --quantMode GeneCounts --outWigType wiggle --outWigStrand Stranded 
#  cmd = " --runMode alignReads --runThreadN 6 --genomeDir {genomeDir} --genomeLoad LoadAndKeep 
--readFilesIn {fq} --readMatesLengthsIn Equal --readFilesCommand 'zcat -1' --clip5pNbases 14 
--outFileNamePrefix {prefix} --outSAMtype BAM Unsorted --outSAMstrandField intronMotif --outSAMreadID Number 
--outFilterMultimapNmax 1 --outFilterMismatchNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 
--outFilterScoreMinOverLread 0.3 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --alignMatesGapMax 200000 --alignIntronMax 200000 
--alignEndsType EndToEnd".format( star=STAR,genomeDir=STAR_INDEX, fq=" ".join(fqs), prefix=prefix )


	"""
	# index=changename_index()
	outfile=pathin+os.path.split(file)[1].split(".")[0]+"/"+os.path.split(file)[1].split(".")[0]
	mkdir(pathin+os.path.split(file)[1].split(".")[0])
	#cmd= "nohup {P1}  --runMode alignReads --runThreadN 6 --genomeDir {P2}  --genomeLoad LoadAndKeep --readFilesIn {P3} --outSAMtype BAM Unsorted --outSAMreadID Number --outFilterMultimapNmax 1 --outFilterMismatchNmax 3 --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3   --outFileNamePrefix {P5} &".format(
	cmd= "{P1}  --runMode alignReads --genomeLoad NoSharedMemory --readFilesCommand 'zcat -1' --clip3pNbases 4 --clip5pNbases 2 --outSAMstrandField intronMotif --runThreadN 6 --genomeDir {P2}    --readFilesIn {P3}  --outFileNamePrefix {P5}  --outSAMreadID Number --outSAMtype BAM SortedByCoordinate --outWigType wiggle --outWigStrand Unstranded --outFilterMultimapNmax 1 --outFilterMismatchNmax 10  --outFilterMismatchNoverLmax 0.1 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 ".format(

		P1=STAR,
		P2=STAR_INDEX,
		P3=file.replace("_&_"," "),
		P4=GTF,
		P5=outfile
		)
	print cmd
	os.system(cmd)

def cufflink(bam):
	sample=os.path.split(bam)[0].split("/")[len(os.path.split(bam)[0].split("/"))-1]
	# print sample

	cmd="{P1}  -G {P2}  --library-type fr-firststrand  --multi-read-correct  --seed 123 -p 12 -q -o {P3} {P4}".format(
	
		P1=CUFFLINKS,
		P2=GTF,
		P3=outpath+count+sample,
		P4=bam
		
		)
	print cmd
	os.system(cmd)

def cuffdiff(bams):
	#cuffdiff --seed 123 --library-type fr-firststrand --FDR 0.05 -b $FA -o O_Y -p 10 $CGTF ../3.Mapping/OA/OA_Aligned.out.bam,../3.Mapping/OB/O
	#B_Aligned.out.bam,../3.Mapping/OC/OC_Aligned.out.bam ../3.Mapping/YA/YA_Aligned.out.bam,../3.Mapping/YB/YB_Aligned.out.bam,../3.Mapping/YC/YC_Aligned.out.bam
	
	cmd="{P1} --seed 123 -b {P5} --library-type fr-firststrand --FDR 0.05 -o {P2} -p 10 {P6} {P3} {P4}".format(
		P1=CUFFDIFF,
		P2=outpath+DEG+"VC_N2",
		P3=bams[0],#control
		P4=bams[1], #treatment
		P5=FA,
		P6=GTF
		)
	print cmd
	# os.system(cmd)

def bam2bed2bg(file):
	# sample=os.path.split(file)[0].split("/")[len(os.path.split(file)[0].split("/"))-1]


	# cmd1="{P1} sort -@ 2 {P2} {P3} | {P1} index {P4}".format(
	# #cmd1="nohup {P1} index {P4} &".format(	
	# 	P1=SAMTOOLS,
	# 	P2=file,
	# 	P3=fileoutpath+sample+".sorted",
	# 	P4=fileoutpath+sample+".sorted.bam",
	# 	)
	# # print cmd1
	# # os.system(cmd1)
	# """
	# infer experiment type
	# infer_experiment.py -r hg19.refseq.bed12 -i Pairend_nonStrandSpecific_36mer_Human_hg19.bam
	# infer_experiment.py -r /picb/molsysbio/usr/caidh/projects/007_bh/data/gencode_v21_lncRNApred_to_R_hg19.bed -i /picb/molsysbio/usr/caidh/projects/013_X/data/1_2014_PLOSONE_PBMC_GSE60216/4_stat_bed_bigwig/SRR1539207.fastq.sorted.bam
	# """
	# cmd3="{P1} -i {P2} -s {P3} -o {P4} --skip-multi-hits --wigsum=1000000000 ".format(#revised the type
	# 	P1=BAM2WIG,
	# 	P2=fileoutpath+sample+".sorted.bam",
	# 	P3=CE10SIZE,
	# 	P4=fileoutpath+sample
	# 	)
	# print cmd3
	# os.system(cmd3)
	cmd4="{P1} -clip {P2} {P3} {P4}".format(
		P1=WIG2BIGWIG,
		P2=file,
		P3=CE10SIZE,
		P4=file.replace(".wig",".bigwig")
		)
	print cmd4
	os.system(cmd4)



def pipeline():
	outpath="/picb/molsysbio2/caidh/032_ASL/GSE101964_celegans_glp1_riboM/"
	# mkdir(outpath+"1_fastq/")
	# sras=glob.glob(outpath+"0_sra/*sra")
	# print sras
	# Parallel(n_jobs=10)(delayed(sra2fastq)(sra,outpath+"1_fastq/") for sra in sras)	

	# fqs=glob.glob(outpath+"1_fastq/*/*fastq.gz")	
	# print fqs
	# Parallel(n_jobs=8)(delayed(STARmappingsingle)(fq,outpath+"2_mapping/") for fq in fqs)	

	wigs=glob.glob(outpath+"2_mapping/*/*wig")
	Parallel(n_jobs=len(wigs))(delayed(bam2bed2bg)(wig) for wig in wigs)



def main(  ):

	pipeline()


if __name__=="__main__":
	main()