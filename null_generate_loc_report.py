import argparse
import RNA
from Bio.Seq import Seq
import shutil
import probability_eval
import json
import os


def nascent_polypeptide(seq, proline_threshold):
	# Checks upstream bases for nascent polypeptide arrest sequence
	# Input: upstream bases, % proline threshold
	# Output: dictionary with info about nascent chain

	dna = Seq(seq)
	amino_acids = dna.translate()
	nascent_info = {}
	nascent_info["Arrest Sequence"] = str(amino_acids)
	return nascent_info

def SD(sequence, candidate_SD_region_len, slippery_loc):
	# Finds the minimum MFE SD interaction in candidate SD region
	# Input: upstream sequeence, anti-SD sequence, candidate region length
	# Output: dictionary with most likely SD interaction

	candidate_SD_region = sequence
	SD_interaction = {}
	for i in range(len(candidate_SD_region) - 4):
		seq = candidate_SD_region[i:i+4]
		if (seq == 'AGGA') or (seq == 'GAGG') or (seq == 'GGAG'):
			SD_interaction = {"SD-like Sequence": seq, "SD start (i)": (slippery_loc - len(candidate_SD_region) + i + 1), "SD end (j)": (slippery_loc - len(candidate_SD_region) + i + 4)}
	return SD_interaction

def find_structure(ref):
	# Looks for downstream structure
	# Input: CaCoFold fold output, original sequence, number of unpaired bases
	# Output: Whether there is a downstream structure that will stall translation, and CaCoFold structure

	# Read in reference
	with open(ref) as f:
		reference = f.read()
	f.close()
	# remove header
	reference = ''.join(reference.split('\n')[1:])
	# compute minimum free energy (MFE) and corresponding structure
	(structure, mfe) = RNA.fold(reference)
	return structure

def structure_info(structure_string):
	# find space and trim string
	space = 0
	for c in range(len(structure_string)):
		if structure_string[c] != '.':
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
		if c == '(':
			p1 += 1
		if c == ')':
			p2 += 1
		ss += c
	paired = p1 + p2
	return len(ss), paired, ss, space

def generate_output(ref, num_bases, upstream_bases, len_SD_region, proline_threshold, report_file, loc, struct_search_start, freqs):
	# Generates dictionary with information about this slippery sequence
	# Only writes to output if there is a downstream structure that will stall translation
	
	structure = find_structure(ref)
	total_len, num_paired, structure_string, space = structure_info(structure)
	slippery_loc = int(loc.split('-')[2])
	with open(upstream_bases) as f:
		reference = f.read()
	f.close()
	upstream= reference[0:len_SD_region]
	# Look for SD-like sequence
	SD_dict = SD(upstream, len_SD_region, slippery_loc)
	# Look for nascent polypeptide arrest sequence
	nascent = nascent_polypeptide(upstream, proline_threshold)
	loc_dict = nascent.copy()
	for key, value in SD_dict.items():
		loc_dict[key] = value
	loc_dict["Slippery Sequence Start (k)"]	= slippery_loc
	loc_dict["Slippery Sequence End (l)"]= slippery_loc + 6
	loc_dict["cDNA"] = loc.split('/')[1]
	loc_dict["Structure"] = structure_string
	loc_dict["h"] = slippery_loc - len_SD_region + 1
	loc_dict["m"] = loc_dict["Slippery Sequence End (l)"] + 4 + space
	loc_dict["o"] = loc_dict["m"] + total_len
	loc_dict["p"] = num_paired
	loc_dict["Full PRF sequence"] = reference[0:25 + space + len(structure_string)]
	with open(ref + ".fullseq.fa", 'w') as f:
		f.write(">Header")
		f.write(loc_dict["Full PRF sequence"])
	f.close()
	loc_dict["LOD"] = probability_eval.calc_probability(loc_dict, freqs)
	with open(report_file, 'a') as f:
		f.write(str(loc_dict))
		f.write("\n")
	f.close()
	return


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', default=None, type=str, help='Path to reference sequence')
	parser.add_argument('-B', default=20, type=int, help='# downstream bases to evaluate for structure')
	parser.add_argument('-U', default=None, type=str, help='Path to upstream bases')
	parser.add_argument('-R', default=20, type=int, help='Length of candidate SD region')
	parser.add_argument('-P', default=None, type=int, help='Proline threshold')
	parser.add_argument('-O', default=None, type=str, help='Path to report file')
	parser.add_argument('-L', default=None, type=str, help='Slippery sequence loc')
	parser.add_argument('-M', default=None, type=str, help='# bases from slippery where structure could begin')
	parser.add_argument('-Q', default=None, type=str, help='path to frequency file')
	args = parser.parse_args()
	generate_output(args.F, args.B, args.U, args.R, args.P, args.O, args.L, args.M, args.Q)

