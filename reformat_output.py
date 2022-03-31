import argparse
import pandas as pd
import json
import matplotlib.pyplot as plt

def num_arresting_residues(string):
	# count number of arresting residues
	n = 0
	for c in string:
		if c == 'P' or c == 'K' or c == 'R':
			n += 1
	return n

def reformat_output(output):
	# Convert report into TSV
	# Input: file with list of dictionaries
	# Output: TSV file with data about each potential PRF site

	dict_list = []
	with open(output, 'r') as f:
		for line in f:
			line = line.replace("'", '"')
			line_dict = json.loads(line)
			dict_list.append(line_dict)
	f.close()
	df = pd.DataFrame(dict_list)
	#move LOD and gene to first and second columns
	first_column = df.pop('gene')
	second_column = df.pop('LOD')
	df.insert(0, 'LOD', second_column)
	df.insert(0, 'gene', first_column)
	df = df.sort_values(by =['LOD'], ascending = False)

	df['Arresting Residues'] = df["Arrest Sequence"].apply(num_arresting_residues)

	df.to_csv(f'{output}.tsv', sep='\t', index=False)
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-O', default=None, type=str, help='Path to report file')
	args = parser.parse_args()
	reformat_output(args.O)
