import argparse
import os

def find_structure(sto):
	# Looks for downstream structure
	# Input: CaCoFold fold output, original sequence, number of unpaired bases
	# Output: Whether there is a downstream structure that will stall translation, and CaCoFold structure
	
	# Read in CaCoFold output
	with open(sto) as f:
		CaCoFold_out = f.read()
	f.close()
	#evaluate structure
	for line in CaCoFold_out.split('#'):
		if line[0:11] == "=GC SS_cons":
			for string in line.split(" "):
				if len(string) > 7:
					structure = string
	return structure

def structure_info(structure_string):
	# find space and trim string
	space = 0
	for c in range(len(structure_string)):
		if structure_string[c] != ':':
			structure_string = structure_string[c:]
			break
		else:
			space += 1
			if space == len(structure_string):
				return 0, 0, structure_string, len(structure_string)

	# find size of structure and num paired
	p1 = 0
	p2 = 0
	ss = str()
	for c in structure_string:
		if p1 == p2 and (p1 != 0):
			break
		if c == '<' or  c == '(':
			p1 += 1
		if c == '>' or  c == ')':
			p2 += 1
		ss += c
	paired = p1 + p2
	return len(ss), paired, ss, space

def generate_output(ref, sto, upstream_bases, loc, len_SD_region):
	# Generates dictionary with information about this slippery sequence
	# Only writes to output if there is a downstream structure that will stall translation
	with open(ref) as f:
		lines = f.read()
	f.close()
	header = lines.split('\n')[0]
	structure = find_structure(sto)
	total_len, num_paired, structure_string, space = structure_info(structure)
	slippery_loc = int(loc.split('-')[2])
	with open(upstream_bases) as f:
		reference = f.read()
	f.close()
	k = slippery_loc
	l = slippery_loc + 6
	h = slippery_loc - len_SD_region + 1
	m = l + 4 + space + 1
	o = m + total_len - 1
	full_sequence = reference[0:o + 1 - h]
	with open(ref + ".fullseq.fa", 'w') as f:
		f.write(header)
		f.write('\n')
		f.write(full_sequence)
	f.close()
	return


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', default=None, type=str, help='Path to reference sequence')
	parser.add_argument('-S', default=None, type=str, help='Path to CaCoFold output')
	parser.add_argument('-U', default=None, type=str, help='Path to upstream bases')
	parser.add_argument('-L', default=None, type=str, help='Slippery sequence loc')
	parser.add_argument('-R', default=20, type=int, help='Length of candidate SD region')
	args = parser.parse_args()
	generate_output(args.F, args.S, args.U, args.L, args.R)

