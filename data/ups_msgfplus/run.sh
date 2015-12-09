#!/bin/bash

msgfplus_bin=../../software/msgfplus/MSGFPlus.jar
dbfile=../fasta/ups.with_contaminants.non_redundant.with_mimic.fasta
dbfile_reversed=../fasta/ups.with_contaminants.non_redundant.with_mimic.reversed.fasta
filename=20080315_CPTAC6_22_6QC1
ms2file=../ups_ms2/${filename}.ms2

java -Xmx3500M -jar ${msgfplus_bin} -s ${ms2file} -d ${dbfile} -o ${filename}.mzid -inst 1 -addFeatures 1
java -Xmx3500M -jar ${msgfplus_bin} -s ${ms2file} -d ${dbfile_reversed} -o ${filename}.reversed.mzid -inst 1 -addFeatures 1
