from Bio import SeqIO
import re

charset=[]
nameset=[]
#whole=[]
with open ("/Users/Ian/Downloads/infile/Feng_infile.fasta", "rU") as f:
    for line in f:
        line=line.strip () # strip the space from before "CHARSET..."
        if "CHARSET" in line:
            line = line.split (" = ") # split the line around the "=", but make careful note of the spaces!
            #locus_name = line.split (" L")
            charset.append (line [1]) # append just the first portion of the split line
            nameset.append (line [0])
            #whole.append (line)
    #print charset
    #print nameset
    #print whole

for c in charset:
    with open (c+"_aligned.fas","w") as out:
        for record in SeqIO.parse ("/Users/Ian/Downloads/infile/Feng_infile.fasta", "fasta"):
            m=re.search("(.*) - (.*);",c) # designate the two portions of the kept line
            if m:
                start = int(m.groups()[0])
                end = int(m.groups()[1])
            #outname = nameset[c]+
            print >> out, ">"+record.id
            print >> out, record.seq [start:end]
