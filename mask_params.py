import argparse
import pandas as pd
import re

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def find_mask_parameters(fa, sto):
	with open(fa) as f:
		reference = f.read()
		reference = reference.split('\n')[1]
	f.close()
	with open(sto) as f:
		align = f.read()
	f.close()
	lines = align.split('\n')
	for line in lines:
		if not line:
			continue
		if line[0] == '#':
			continue
		else:
			top_hit = re.split(r'\s\s*', line)[1]
			break
	if len(top_hit) > len(reference):
		top_hit = top_hit[:10]
		hamming_distances = []
		for i in range(len(reference) - len(top_hit)):
			ref = reference[i:i+len(top_hit)]
			hamming_distances.append((i, hamming_distance(top_hit, ref)))
		min_ = hamming_distances[0]
		for tup in hamming_distances[1:]:
			if tup[1] < min_[1]:
				min_ = tup
		print(13 - min_[0])
		return
	elif len(reference) == len(top_hit):
		print(13)
		return
	else:
		hamming_distances = []
		for i in range(len(reference) - len(top_hit)):
			ref = reference[i:i+len(top_hit)]
			hamming_distances.append((i, hamming_distance(top_hit, ref)))
		min_ = hamming_distances[0]
		for tup in hamming_distances[1:]:
			if tup[1] < min_[1]:
				min_ = tup
		print(13 - min_[0])
		return

if __name__ == "__main__":
 	parser = argparse.ArgumentParser()
 	parser.add_argument('-F', default=None, type=str, help='Original sequence')
 	parser.add_argument('-S', default=None, type=str, help='Alignment')
 	args = parser.parse_args()
 	find_mask_parameters(args.F, args.S)


