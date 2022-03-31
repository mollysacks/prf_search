import argparse
import pandas as pd
import json
from Bio.Seq import Seq
import re
import numpy as np
import math
import RNA

def to_RNA(dna):
	rna = ""
	for c in dna:
		if c == 'T':
			rna = rna + 'U'
		else:
			rna = rna + c
	return rna

def null(seq, freq):
	# Calculate null sequence probability
    prob = 0
    for base in seq:
        prob += np.log2(freq[base])
    return 2 ** (prob)

def calc_upstream_prob(seq, freq, tildes, stars, arresting):
	# calculate upstream probability
    F = upstream_subset(seq[0:15])
    prob = stars[F] * null(seq[0:15], freq) / tildes[F]
    return(prob)

def upstream_subset(seq):
	# parses upstream sequence for features
	arresting = ['K', 'R', 'P']
	if len(seq) != 15:
		seq = seq[0:15]
	np = Seq(seq).translate()
	arresting_bases = 0
	for a in np:
		if a in arresting:
			arresting_bases += 1
	if ('AGGA' in seq) or ('GGAG' in seq) or ('GAGG' in seq): 
		if arresting_bases >= 2:
			return(f'SD_NP{arresting_bases}')
		else:
			return('SD')
	elif arresting_bases >= 2:
		return(f'NP{arresting_bases}')
	else:
		return('Neither')

def slippery_prob(seq, freq, ss):
	# calculates slippery sequence probability
	if seq[0] != seq[1]:
		return null(seq[1:], freq) * (2 ** -6)
	else:
		return (null(seq[1:], freq)/ ss[seq[0]])

def P_DS(sequence, freq, structure_priors, l, n, structure):
	# Calculates probability of observing an RNA secondary structure taht could inhibit the ribosome in the downstream sequence
	m_dist = list()
	if len(structure) <= 10:
		return(null(sequence[l+1:], freq) * (2 ** (-(len(sequence) - l + 1))))
	min_marg = 1
	for value in structure_priors.values():
		if int(value) == 0:
			continue
		if int(value) < min_marg:
			min_marg = int(value)
	for m in range(l+6,len(sequence)): #start at 6th base after end of slippery sequence
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
		# if there are no structures, continue to next value of m
		if len(subopts) == 0:
			continue
		for s in subopts:
			if s.energy <= 0:
				total_mfe += s.energy
				if s.structure[0] == '(': # check if there is structure that starts at base m
					struct_mfe += s.energy
		# if there are no MFE strctures, continue to next value of m
		if total_mfe == 0:
			continue
		if struct_mfe == 0:
			continue
		likelihood = struct_mfe/total_mfe
		prior = null(sequence[l+1:], freq)
		marg = structure_priors[str(len(seq))]
		if marg == 0:
			# use minimum marginalization to avoid div/0
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

def calc_probability(site_info, freqs, reference):
	h = site_info['h']
	# import other loci, change to 0 indexing
	if 'SD start (i)' in site_info:
		i = site_info['SD start (i)'] - h
		j = site_info['SD end (j)'] - h
	else:
		i = 0
		j = 0
	k = site_info['Slippery Sequence Start (k)'] - h 
	l = site_info['Slippery Sequence End (l)'] - h
	m = site_info['m'] - h 
	n = site_info['n'] - h 
	#set h to zero
	h = h - h
	
	structure = site_info['Structure']
	seq = reference
	# calculate index probabilities
	# load frequency data
	with open(freqs, 'r') as f:
		info = f.read()
		priors = info.replace("'", "\"")
		priors = json.loads(priors)
	f.close()

	nucleotide_freq = priors['freq']
	structure_marg = priors['Structure params']
	tildes = priors['tildes']
	stars = priors['stars']
	ss = priors['SS_t']

	# calculate probabilities
	upstream_prob = calc_upstream_prob(seq, nucleotide_freq, tildes, stars, ['R', 'P', 'K'])
	slip_prob = slippery_prob(seq[k:l+1], nucleotide_freq, ss)
	struct_prob = P_DS(reference, nucleotide_freq, structure_marg, l, n, site_info['Structure'])	
	lodd = {}
	lodd['Structure LOD'] = np.log2(struct_prob/null(seq[l+1:], nucleotide_freq))
	lodd['Upstream LOD'] = np.log2(upstream_prob/null(seq[h:k+1], nucleotide_freq))
	lodd['Slippery LOD'] = np.log2(slip_prob/null(seq[k+1:l+1], nucleotide_freq))
	total_prob = upstream_prob * slip_prob * struct_prob

	null_prob = null(seq, nucleotide_freq)
	#lodd['subopt_structs'] = subopt_structs
	lodd['LOD'] = np.log2(total_prob/null_prob)
	#lodd['Structure MFE'] = mfe
	return lodd

