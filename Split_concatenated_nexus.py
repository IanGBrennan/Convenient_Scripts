from Bio.Nexus import Nexus

aln = Nexus.Nexus()
#aln.read('my-properly-formatted-nexus-file.nex')
aln.read('/Users/Ian/Desktop/G.japonicus.RELEC.nex')

# assuming your partitions are defined in a charset block like:
#
# begin sets;
# charset bag2 = 1-186;
# charset bag3 = 187-483; 
# charset bche = 484-990;
# end;

# get count of charsets:
len(aln.charsets.keys())

# take a gander at the charsets:
aln.charsets()

# split the concatenated file into charsets files (prepending whatever text you place after filename='')
aln.write_nexus_data_partitions(filename='my', charpartition=aln.charsets)

# this will output in the os.getcwd():
# 
# test_bag2
# test_bag3
# test_bche
#os.getcwd()