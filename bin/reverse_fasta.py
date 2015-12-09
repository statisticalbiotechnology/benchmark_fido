#!/usr/bin/python

def readFasta(fp):
  name, seq = None, []
  for line in fp:
    line = line.rstrip()
    if line.startswith(">"):
      if name: yield (name, "".join(seq))
      name, seq = line[1:].split(" ")[0], []
    else: seq.append(line)
  if name: yield (name, "".join(seq))

def reverseFasta(inputFile, outputFile):
  with open(inputFile, 'rb') as f:
    with open(outputFile, 'wb') as w:
      for proteinName, proteinSeq in readFasta(f):
        print >> w, ">random_" + proteinName
        print >> w, proteinSeq[::-1]
        
inputFile = "../data/fasta/ups.with_contaminants.non_redundant.with_mimic.fasta"
outputFile = "../data/fasta/ups.with_contaminants.non_redundant.with_mimic.reversed.fasta"

reverseFasta(inputFile, outputFile)
