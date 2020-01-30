#! /usr/bin/env python
import sys
import re
''' Description: takes phylip files (given by user) either sequencial or interleaved 
and checks for empty sequences (i.e. not containing any ATCG) and deletes them. The 
reduced dataset is written in the directory as trimmed.phy
'''
print(sys.argv[1])
output_name=[]
input_name=sys.argv[1]
split_up=input_name.split(".")
out_name=split_up[0]


f1 = open(sys.argv[1], 'r')
f2 = open(out_name+'_trimmed.phy', 'w')
content = f1.readlines()
new = []
pattern = re.compile(r'[ATCG]')
for elem in range(len(content)):
	new.append(content[elem].split())
for elem1 in range(1, len(new)):
	if len(new[elem1]) == 2:
		if pattern.findall(new[elem1][1]):
			pass
		else:
			new[0][0] = int(new[0][0]) - 1
f2.write(str(new[0][0]) + " ")
f2.write(new[0][1] + "\n")
for elem1 in range(1, len(new)):
	if len(new[elem1]) == 2:
		if pattern.findall(new[elem1][1]):
			f2.write(new[elem1][0] + " ")
			f2.write(new[elem1][1] + "\n")						
f2.write("\n")			
for elem1 in range(len(new)):
	if len(new[elem1]) == 1:
		if pattern.findall(new[elem1][0]):
			f2.write(new[elem1][0] + "\n")		
f1.close()
f2.close()
