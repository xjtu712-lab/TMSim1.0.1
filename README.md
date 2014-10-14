
Manual of TMSim 1.0.1

Contact:Yu Geng, gengyu@std.xjtu.edu.cn
		Zhongmeng Zhao, zmzhao@mail.xjtu.edu.cn
		2014,March


Introduction

TMSim is a reads simulator for the multi-level structure of the cancer clone.This program is designed to embody the characteristic of interdependent and interaction between different level subclones that can emerge inheritance and variability.Currently researches usually adopt the joint analysis of tumor-normal pair data in cancer.Hence,this program can generate *.sim files and set variations of sites about normal and subclones.It can random sampling to reads based on *.sim files of normal and subclones and form the simulation environment including all kinds of variations and generate two files of left.fq and right.fq that exit in positive and reverse chain,respectively. 

Installation

1.TMSim aims for human diploid genomes. It runs on 64-bit Linux system .
2.You can download and unpack using "tar -xzvf  <destination folder>.tar.gz" and execute directly object file.
3.Or download the source code, unpack to ${destination folder} with the method above, and compile by using "g++ <filename.cpp>" at ${destination folder}/TMSim.Then execute directly.


How to use TMSim program package

1.NorSim:the Normal Simulator for the data generator(program: norsim.cpp)

NorSim set variations of site for nomal cell sequence.Run this program must have provided two file names including <refrence.fa>(fasta format,human reference genome) and <nor.sim>.

Usage:   norsim [options] <ref.fa> <nor.sim>
Options: -r FLOAT		rate of mutation [0.0010000000]
         -R FLOAT      		fraction of indels [  0.100000]
         -X FLOAT      		probability an indel is extended [  0.300000]
         -D FLOAT      		delete rate in indel [  0.500000]
         -B FLOAT      		BB rate in mutation [  0.333330]
         -I <delpos.txt>  		input long indel set file
         -A <nor_AB.idx>  		output the positions of AB mutation in normal
         -o <nor_simout.txt>  	result for runing case

For example: norsim  -r 0.01 -R 0.00 -X 0.00  -A nor_AB.idx  chr1.fa nor.sim

Notes:  an option"-A <nor_AB.idx>" that records the positions of AB mutation gype in normal, if  continue to generate *.sim of subclones based on normal,we should select this option.Otherwise,you may not choose it .

2.TumSim:the Subclones Simulator for the data generator(program: tumsim.cpp)

Once nor.sim and nor_AB.idx files are available, a way to generate the *.sim file of subclones in tumor is:
Usage:   tumsim [options] <ref.fa> <base.sim> <nor_AB.idx> <subclone.sim>
Options: -r FLOAT      		rate of mutation [0.0000240000]
         -R FLOAT      		fraction of indels [  0.100000]
         -X FLOAT      		probability an indel is extended [  0.300000]
         -D FLOAT      		delete rate in indel [  0.500000]
         -B FLOAT      		BB rate in mutation [  0.333330]
         -A FLOAT      		normal AB rate in mutation [  0.500000]
         -p FLOAT      		the position of high desity [0.500000]
         -n int      		mutation number in high desity [120]
         -l int       		length of high desity [0]
         -I <delpos.txt>  		input long indel set file
         -N <other_chged.idx> 	input changed positions in other subclone(may multi-choices)
         -C <*_chged.idx>  	output changed positions in this turmor
         -o <*_simout.txt>  	output result for runing case

For example: 
	(1)tumsim [options]  chr1.fa nor.sim nor_AB.idx s1.sim
	(2)tumsim [options] -C s3_chged.idx chr1.fa s1.sim nor_AB.idx s3.sim
    (3)tumsim [options] -N s3_chged.idx chr1.fa s1.sim nor_AB.idx s2.sim
Notes:
If s1.sim only inherites variations of normal and generate itself variation sits,then use the first command form.  
If s1 produces separately two subclones that are s2 and s3, the variations of s3 should be avoided by s2 ,you should first produce s3.sim including an option parameter of "-C s3_chged.idx",then produce s2.sim including an option parameter of "-N s3_chged.idx".

3.ReadGen:the Reads Generator from *.sim (program: readgen.cpp)

when all the *.sim files are obtained,random generate two files of left.fq and right.fq that exit in positive and reverse chain,respectively.
Usage:   readgen [options] <ref.fa> <*.sim> <left.fq> [<right.fq>]
Options: -d INT		outer distance between the two ends [10000]
         -s INT        	standard deviation [10]
         -l INT        	length of left  read [100]
         -r INT        	length of right read [100]
         -c float      	cover rate of pair_ends read [0.500000]
         -e float      	error rate in generate paired_end reads [0.010000]
         -S       		output single read(not paired_end reads)
         -k       		keep 'N' character in reads
         -I <delpos.txt>  	long indel set file
         -O <genread.log>  	result for runing case

For example: readgen  chr1.fa nor.sim nor_left.fq nor_right.fq
            
Notes:
1.Assuming diploid individuals;
2.This programe supports single reads and paired-end reads that can be achieved whether or not the parameter"-S" are selected ,respectively.
					
				



