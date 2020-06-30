#!/bin/bash

module add UHTS/Analysis/seqtk/1.2

seqtk sample -s100 E1_S1_L001_R1_001.fastq.gz 20000 > E1_S1_L001_R1_001-SUB.fastq.gz
seqtk sample -s100 E1_S1_L001_R2_001.fastq.gz 20000 > E1_S1_L001_R2_001-SUB.fastq.gz
seqtk sample -s100 E1_S1_L002_R1_001.fastq.gz 20000 > E1_S1_L002_R1_001-SUB.fastq.gz
seqtk sample -s100 E1_S1_L002_R2_001.fastq.gz 20000 > E1_S1_L002_R2_001-SUB.fastq.gz

seqtk sample -s100 A3_S2_L001_R1_001.fastq.gz 20000 > A3_S2_L001_R1_001-SUB.fastq.gz
seqtk sample -s100 A3_S2_L001_R2_001.fastq.gz 20000 > A3_S2_L001_R2_001-SUB.fastq.gz
seqtk sample -s100 A3_S2_L002_R1_001.fastq.gz 20000 > A3_S2_L002_R1_001-SUB.fastq.gz
seqtk sample -s100 A3_S2_L002_R2_001.fastq.gz 20000 > A3_S2_L002_R2_001-SUB.fastq.gz