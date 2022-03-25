import argparse
import matplotlib.pyplot as plt
import re
import pandas as pd

def find_distribution(fasta):
	with open(fasta, 'r') as f:
		lines = f.read()
	f.close()
	genes_list = lines.split('\n>')
	lengths = []
	for gene in genes_list:
		sequence = re.sub(r'^.*?\n', '', gene)
		lengths.append(len(sequence))
	len_df = pd.DataFrame(lengths)
	mean = len_df.mean()
	var = len_df.var()
	return mean, var

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', default=None, type=str, help='Path to input cDNA fasta')
	args = parser.parse_args()
	find_distribution(args.F)
