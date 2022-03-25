import argparse
import pandas as pd
import json
import matplotlib.pyplot as plt
import csv
from probability_eval import upstream_subset
import matplotlib.lines as mlines

def num_arresting_residues(string):
	n = 0
	for c in string:
		if c == 'P' or c == 'K' or c == 'R':
			n += 1
	return n

def colors(df):
	cols = list([plt.cm.tab10(i) for i in range(10)]) 
	categories = ['Neither', 'SD', 'NP2', 'NP3', 'NP4', 'NP5', 'SD_NP2', 'SD_NP3', 'SD_NP4', 'SD_NP5']
	assignments = {}
	for i in range(10):
		assignments[categories[i]] = cols[i]
	df['upstream'] = df['Full PRF sequence'].apply(upstream_subset)
	df['color'] = df['upstream'].apply(lambda x: assignments.get(x))
	return df, assignments

def merge(out_df, no_cacofold):
	#2163 out
	#2349
	conserved = list(zip(out_df['cDNA'], out_df['Slippery Sequence Start (k)']))
	for i, row in no_cacofold.iterrows(): 
		if (row['gene'], row['Slippery Sequence Start (k)']) not in conserved:
			row['Structure'] = ''
			out_df = out_df.append(row)
	return out_df

def create_plots(output, null, no_cacofold):
	with open(output) as f:
		out_df = pd.read_csv(f, sep="\t")
	f.close()
	with open(null) as f:
		null_df = pd.read_csv(f, sep="\t")
	f.close()
	with open(no_cacofold) as f:
		no_cacofold = pd.read_csv(f, sep="\t")
	out_df = merge(out_df, no_cacofold)

	
	conserved = out_df.loc[out_df['h'] != None]
	out_df = out_df.loc[out_df['Len Vienna region'] == 80]
	out_df['Struct len'] = out_df['Structure'].apply(len)
	
	out_df, assignments = colors(out_df)
	print(len(out_df))
	LOD_mean = out_df['LOD'].mean()
	LOD_median = out_df['LOD'].median()
	LOD_var = out_df['LOD'].var()
	print(LOD_mean, LOD_median, LOD_var)
	out_df['Slippery Sequence'] = out_df['Full PRF sequence'].str.slice(14,21)
	out_df.to_csv(f'{output}.final.tsv', sep='\t', index=False)

	plt.hist(out_df['LOD'], range=[-15,15], bins = 80)
	plt.xlim(-15, 15)
	plt.xlabel('E. coli LOD score')
	plt.ylabel('Number of occurrences')
	plt.title('E. coli LOD score distribution')
	plt.show()

	plt.hist(null_df['LOD'], range=[-15,15], bins = 80)
	plt.xlim(-15, 15)
	plt.xlabel('LOD score')
	plt.ylabel('Number of occurrences')
	plt.title('Null LOD score distribution')
	plt.show()

	artists = []
	cats = []
	for cat, col in assignments.items():
		p = mlines.Line2D([], [], color=col, marker='o', linestyle='None',
                          markersize=10, label=cat)
		artists.append(p)
		cats.append(cat)
	plt.scatter(out_df['Structure LOD'], out_df['LOD'], alpha=.8, c=out_df['color'], label=out_df['color'])
	plt.xlim(-25, 5)
	plt.ylim(-15,15)
	plt.xlabel('Structure LOD')
	plt.ylabel('LOD')
	plt.title('E. coli LOD vs. Structure LOD')
	plt.legend(artists, cats)
	plt.show()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-O', default=None, type=str, help='Path to tsv')
	parser.add_argument('-N', default=None, type=str, help='Path to null tsv')
	parser.add_argument('-C', default=None, type=str, help='Path to no CaCoFold tsv')
	args = parser.parse_args()
	create_plots(args.O, args.N, args.C)
