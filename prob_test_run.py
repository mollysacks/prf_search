import pandas as pd
import json
import probability_eval

dict_list = []
with open('all2022_01_10-11_28_49.report', 'r') as f:
	for line in f:
		line = line.replace("'", '"')
		line_dict = json.loads(line)
		dict_list.append(line_dict)
f.close()
df = pd.DataFrame(dict_list)
#d = {'Arrest Sequence': 'GATKA', 'SD-like Sequence': 'GGAG', 'SD start (i)': 1268, 'SD end (j)': 1271, 'SD MFE': -9.300000190734863, 'Slippery Sequence Start (k)': 1283, 'Slippery Sequence End (l)': 1289, 'cDNA': 'MG1655', 'Structure': '<<<<<<-<<<<<_____>>>->>->>>>>>', 'm': 1295, 'n': 1306, 'o': 1324, 'h': 1269, 'Full PRF sequence': 'GGAGCAACCAAAGCAAAAAAGAGTGAACCGGCAGCCGCTACCCGCGCGCGGCCGGT', 'LODD': 0.006043487666226613}


for d in dict_list:
	probability_eval.calc_probability(d)



# def pairs(i,j):
# 	if i == 'A' and j == 'T':
# 		return (freq['A'] + freq['T'])/2
# 	if i == 'T' and j == 'A':
# 		return (freq['A'] + freq['T'])/2
# 	if i == 'G' and j == 'C':
# 		return (freq['G'] + freq['C'])/2
# 	if i == 'C' and j == 'G':
# 		return (freq['G'] + freq['C'])/2
# 	return 0

# def new_structure_prob(sequence, num_pairs):
# 	i = 0
# 	rem_bases = len(sequence)
# 	rem_paired = 2*num_pairs
# 	rem_unpaired = len(sequence) - rem_paired
# 	prob = 1
# 	while sequence:
# 		i = sequence[0]
# 		j = sequence[-1]
# 		if pairs(i,j) > 0:
# 			prob *= pairs(i,j) * (rem_paired/rem_sequence) 
# 			sequence = sequence[1:-1]
# 			rem_paired -= 2
# 			rem_sequence -= 2
# 		else:


# 	return None



# new_structure_prob("AAATGGTATTCCTGTCCACGATACCAT", 9)