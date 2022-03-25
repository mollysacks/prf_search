import argparse
import numpy.random
import random
import re
import pandas as pd
import numpy as np
import probability_eval
import RNA
import json
from Bio.Seq import Seq

def null(seq, freq):
	# Calculate null sequence probability
    prob = 0
    for base in seq:
        prob += numpy.log2(freq[base])
    return 2 ** (prob)

def generate_lengths(n, fasta):
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
	sd = len_df.std()
	i = 0
	rand_lengths = []
	for i in range(n):
		length = -1
		while length < 1:
			length = int(numpy.random.normal(mean, sd))
		rand_lengths.append(length)
	return rand_lengths

def generate_random_genes(lengths, freq):
	cdnas = {}
	bases = []
	weights = []
	for base, frequency in freq.items():
		bases.append(base)
		weights.append(frequency)
	i = 0
	j = 0
	for l in lengths:
		rand_seq = ''.join(random.choices(bases, weights = weights, k=l))
		cdnas[f'RandomGene{str(i)}'] = rand_seq
		i += 1
	return cdnas

def slippery_sequence_search(genes, gene_lengths):
	length = 91
	pattern_matches = []
	for name, seq in genes.items():
		if name in gene_lengths:
			length = gene_lengths[name]
		elif len(seq) < 4000:
			length = len(seq)
		else:
			length = 4000
		segments = list() #initialize list of segments
		i = 14 # start at last base of fifth codon
		while (i <= length - 80): #len(reference) - 7): # go until last 11
			if re.match(r'([ATCG])\1{2}([AT])\2{2}[ATCG]', seq[i:i+7]): # check if slippery sequence
				if len(seq) - i > 91: # if there are > length bases left
					new_seq = seq[i - 14: i + 92] # include an extra 12 bases in case hmmer cuts some off
				else: # if there are < length bases left
					new_seq = seq[i - 14:] #extra bases will be trimmed later on
				pattern_matches.append({'gene': name, 'seq': new_seq, 'k': i})
			i = i + 3 # go to next codon
	return pattern_matches

def vienna_best_structure(seq):
	SS, MFE = RNA.fold(seq)
	n = 25
	up = False
	down = False
	struct = False
	structure = ''
	for c in SS:
		n += 1
		if c == '(':
			if down:
				n -= 1
				break
			struct = True
		if struct:
			structure += c
		if c == ')':
			down = True
	extras = 0
	for c in reversed(structure):
		if c != '.':
			break
		else:
			extras += 1
	n = n - extras
	structure = structure[:-extras]
	return n, structure

def P_DS(sequence, freq, structure_priors, l, n, structure):
	m_dist = list()
	if len(structure) <= 10:
		return(null(sequence[l+1:], freq) * (2 ** (-(len(sequence) - l + 1))))
	min_marg = 1
	for value in structure_priors.values():
		if int(value) == 0:
			continue
		if int(value) < min_marg:
			min_marg = int(value)
	for m in range(l+6,len(sequence)):
		space = len(sequence[l+1:m])
		if (space > 4) and (space < 10):
		    pm = .1999
		else:
		    s = space - 9
		    pm = .0005 * (2 ** -s)
		seq = sequence[m:]
		total_mfe = 0
		struct_mfe = 0
		if len(seq) == 0:
			continue
		a = RNA.fold_compound(seq)
		subopts = list(a.subopt(100))
		if len(subopts) == 0:
			continue
		for s in subopts:
			if s.energy <= 0:
				total_mfe += s.energy
				if s.structure[0] == '(':
					struct_mfe += s.energy
		if total_mfe == 0:
			continue
		if struct_mfe == 0:
			continue
		likelihood = struct_mfe/total_mfe
		prior = null(sequence[l+1:], freq)
		marg = structure_priors[str(len(seq))]
		if marg == 0:
			m_dist.append((pm, likelihood * prior/ min_marg))
			continue
		posterior = (likelihood*prior)/marg
		m_dist.append((pm, posterior))
	P_DS = 0
	if len(m_dist) == 0:
		return null(seq[l+1:], freq) * (2**(-len(sequence[l+1:])))
	for Pm, DS_m in m_dist:
	    P_DS += Pm * DS_m
	return P_DS

def prob_eval(pattern_matches, nucleotide_freq, tildes, stars, structure_marg, ss):
	records = []
	for match in pattern_matches:
		lodd = {}
		lodd['gene'] = match['gene']
		seq = match['seq']
		lodd['Slippery Sequence Start (k)'] = match['k']
		upstream = seq[0:15]
		dna = Seq(upstream)
		amino_acids = dna.translate()
		lodd["Arrest Sequence"] = str(amino_acids)
		lodd['Full PRF sequence'] = seq
		if ('AGGA' in upstream) or ('GGAG' in upstream) or ('GAGG' in upstream):
			lodd['SD'] = True
		else:
			lodd['SD'] = False 
		upstream_prob = probability_eval.calc_upstream_prob(seq, nucleotide_freq, tildes, stars, ['R', 'P', 'K'])
		slip_prob = probability_eval.slippery_prob(seq[14:21], nucleotide_freq, ss)
		n, struct = vienna_best_structure(seq[26:])
		struct_prob = P_DS(seq, nucleotide_freq, structure_marg, 20, n, struct)	
		lodd['Structure LOD'] = np.log2(struct_prob/null(seq[21:], nucleotide_freq))
		lodd['Upstream LOD'] = np.log2(upstream_prob/null(seq[0:15], nucleotide_freq))
		lodd['Slippery LOD'] = np.log2(slip_prob/null(seq[15:21], nucleotide_freq))
		total_prob = upstream_prob * slip_prob * struct_prob
		lodd['n'] = n
		lodd['len vienna region'] = len(seq) - 26
		null_prob = null(seq, nucleotide_freq)
		#lodd['subopt_structs'] = subopt_structs
		lodd['LOD'] = np.log2(total_prob/null_prob)
		records.append(lodd)
	df = pd.DataFrame(records)
	return df

def random_LODs(n, fasta, inp):
	with open('gene_lengths_e_coli') as f:
		glength = f.read()
		glength = glength.replace("'", '"')
		glength = json.loads(glength)
		f.close()
	with open('frequencies', 'r') as f:
		info = f.read()
		priors = info.replace("'", "\"")
		priors = json.loads(priors)
	f.close()
	freq = priors['freq']
	structure_marg = priors['Structure params']
	tildes = priors['tildes']
	stars = priors['stars']
	ss = priors['SS_t']
	if inp == None:
		lengths = generate_lengths(n, fasta)
		genes = generate_random_genes(lengths, freq)
		out = "null.tsv"
	else:
		out = f'{inp}.results.tsv'
		with open(inp) as f:
			seqs = f.read()
		f.close()
		genes_list = seqs.split('\n>')
		genes = {}
		for gene in genes_list:
			header = gene.split('\n')[0]
			seq = ''.join(gene.split('\n')[1:])
			gene_name = header.split(' ')[0]
			genes[gene_name] = seq
	pattern_matches = slippery_sequence_search(genes, glength)
	df = prob_eval(pattern_matches, freq, tildes, stars, structure_marg, ss)
	df = df.sort_values(by =['LOD'], ascending = False)
	df.to_csv(out, sep='\t', index=False)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-N', default=1000, type=int, help='Number of genes to generate')
	parser.add_argument('-F', default=None, type=str, help='Fasta to imitate gene length distribution')
	parser.add_argument('-I', default=None, type=str, help='Optional specify input file')
	args = parser.parse_args()
	random_LODs(args.N, args.F, args.I)
