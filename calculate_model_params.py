from Bio import Seq
import RNA
import argparse
import random
import re
import numpy

def null(seq, freq):
	# Calculate null sequence probability
    prob = 0
    for base in seq:
        prob += numpy.log2(freq[base])
    return 2 ** (prob)

def all_sequences(L, freq, seq_list):
	# Recursively generate all possible sequences of length L
    if L == 0:
        return seq_list
    else:
        new_seq_list = []
        for seq in seq_list:
            new_seq_list.append(seq + 'A')
            new_seq_list.append(seq + 'C')
            new_seq_list.append(seq + 'G')
            new_seq_list.append(seq + 'T')
        return all_sequences(L-1, freq, new_seq_list)

def calculate_upstream_tildes(L, freq):
	# Calculate expected probability of drawing a sequence with 
	# a given set of features from the null distribution of 
	# sequences of length L
	i = 0
	bases = ["A","G","C","T"]
	weights = [freq['A'], freq['G'], freq['C'], freq['T']]
	arresting = ['K', 'P', 'R']
	max_arresting = int(L / 3)
	totals = {'Neither': 0, 'SD': 0}
	for i in range(max_arresting + 1):
		if i < 2:
			continue
		totals[f'SD_NP{i}'] = 0
		totals[f'NP{i}'] = 0
	total_null = 0
	i = 0
	while i < 100000: # generate 10,000,000 sequences
		random_sequence = ''.join(random.choices(bases, weights = weights, k=15))
		np = Seq.Seq(random_sequence).translate()
		arresting_bases = 0
		total_null += null(random_sequence, freq)
		for a in np:
			if a in arresting:
				arresting_bases += 1
		# Shine-Dalgarno like sequences
		if ('AGGA' in random_sequence) or ('GGAG' in random_sequence) or ('GAGG' in random_sequence): 
			if arresting_bases >= 2:
				totals[f'SD_NP{arresting_bases}'] += null(random_sequence, freq)
			else:
				totals['SD'] += null(random_sequence, freq)
		# Nascent peptide arrest sequence
		elif arresting_bases >= 2:
			totals[f'NP{arresting_bases}'] += null(random_sequence, freq)
		# Neither
		else:
		    totals['Neither'] += null(random_sequence, freq)
		i += 1
	tilde = {}
	for key, value in totals.items():
		tilde[key] = totals[key]/total_null
	return(tilde)

def calculate_stars(bits, tildes):
	# Calculate parameters with which to score each upstream feature
    denom = 1
    stars = {}
    for feature, feature_tilde in tildes.items():
        if feature == 'Neither':
            continue
        denom += ((2 ** (bits[feature])) * feature_tilde)/tildes['Neither']
    neither_star = 1/ denom
    for feature, feature_tilde in tildes.items():
        if feature == 'Neither':
            stars[feature] = neither_star
        else:
            stars[feature] = (2 ** bits[feature]) * neither_star * feature_tilde / tildes['Neither']
    return stars

def calc_structure_prior(freq):
	# Calculate expected probability of observing an RNA
	# structure that starts at the first base in the sequence
	# For sequences of length 10-76
	print("Calculating Structure Prior")
	dna = ["A","G","C","T"]
	weights = [freq['A'], freq['G'], freq['C'], freq['T']]
	avg_prob = 0
	structure_priors = {}
	for k in range(1, 85):
		print(k)
		i = 0
		while i <= 2000:
			random_sequence = ''.join(random.choices(dna, weights = weights, k=k))
			a = RNA.fold_compound(random_sequence)
			subopts = list(a.subopt(100))
			total_energy = 0
			m_start_energy = 0
			for s in subopts:
				if s.energy <= 0:
					total_energy += s.energy
					if s.structure[0] == '(':
						m_start_energy += s.energy
			if total_energy == 0:
				avg_prob = 0
			else:
				avg_prob += m_start_energy/total_energy
			i += 1
		avg_prob = avg_prob/1000
		structure_priors[str(k)] = avg_prob
	print("Done calculating structure prior")
	return structure_priors

def codon_frequencies(fa):	 	
	with open(fa) as f:
		combined_fasta = f.read()
	f.close()
	cDNAs_list = combined_fasta.split('\n>')
	cDNAs_raw_list = []
	for cDNA in cDNAs_list:
		raw_cDNA = ""
		lines = cDNA.split('\n')
		if len(lines) == 1:
			continue
		for line in lines[1:]:
			raw_cDNA = raw_cDNA + line
		cDNAs_raw_list.append(raw_cDNA)
	nucleotides = {'A':0, 'T': 0, 'G': 0, 'C': 0}
	for cDNA in cDNAs_raw_list:
		if cDNA == 'NIL':
			continue
		for char in cDNA:
			nucleotides[char] = nucleotides[char] + 1
	nucleotide_frequencies = {}
	for nucleotide, count in nucleotides.items():
		nucleotide_frequencies[nucleotide] = nucleotides[nucleotide] / sum(nucleotides.values())
	structure_prior = calc_structure_prior(nucleotide_frequencies)
	tildes = calculate_upstream_tildes(15, nucleotide_frequencies)
	bit_score = {}
	bit_score['Neither'] = 0
	bit_score['SD'] = 2
	bit_score['NP2'] = 1
	bit_score['NP3'] = 2
	bit_score['NP4'] = 3
	bit_score['NP5'] = 4
	bit_score['SD_NP2'] = 3
	bit_score['SD_NP3'] = 4
	bit_score['SD_NP4'] = 5
	bit_score['SD_NP5'] = 6
	stars = calculate_stars(bit_score, tildes)
	params = {}
	SS_t = {'A':0, 'C':0, 'G':0, 'T':0}
	hexamers = all_sequences(6, nucleotide_frequencies, [''])
	for seq in hexamers:
	    if re.match(r'([ATCG])\1{1}([AT])\2{2}[ATCG]', seq):
	        SS_t[seq[0]] += null(seq, nucleotide_frequencies)
	params['SS_t'] = SS_t
	params['Structure params'] = structure_prior
	params['freq'] = nucleotide_frequencies
	params['tildes'] = tildes
	params['stars'] = stars
	return params

def write_output(params):
	with open('frequencies', 'w') as f:
		f.write(str(params))
	f.close()
	return

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', default=None, type=str, help='Path to fasta')
	args = parser.parse_args()
	params = codon_frequencies(args.F)
	write_output(params)
