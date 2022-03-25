import argparse
import pandas as pd
import re

def extract_conservation_info(loc_dict, cinfo):
	# read in conservation info
	with open(cinfo) as f:
		c_info = f.read()
	f.close()
	records = {}
	for line in c_info.split("\n")[:-1]:
		if (line[0] != '#') and (line[0] != '/'):
			record = {}
			info = re.split(r'\s\s*', line)
			total = sum([float(i) for i in info[2:]])
			record['A'] = float(info[2])
			record['C'] = float(info[3])
			record['G'] = float(info[4])
			record['T'] = float(info[5])
			records[int(info[1]) - 1] = record # for 0 indexing]

	# import sequence and feature info
	h = loc_dict['h']
	# import other loci, change to 0 indexing
	if 'SD start (i)' in loc_dict.keys():
		i = loc_dict['SD start (i)'] - h
		j = loc_dict['SD end (j)'] - h 
	else:
		i = 1234
		j = 1234
	k = loc_dict['Slippery Sequence Start (k)'] - h
	l = loc_dict['Slippery Sequence End (l)'] - h
	m = loc_dict['m'] - h
	n = loc_dict['n'] - h
	#set h to zero
	h = 0
	seq = loc_dict['Full PRF sequence']
	up_cons = 0
	ss_cons = 0
	ds_cons = 0
	total_cons = 0
	sd_cons = 0
	a = 0
	while a <= k:
		if (a >= i) and (a <= j):
			sd_cons += records[a][seq[a]]
		up_cons += records[a][seq[a]] 
		total_cons += records[a][seq[a]]
		a += 1
	up_cons /= 15
	sd_cons /= 4

	# #go back one to check slippery sequence
	a -= 1
	# subtract last iteration
	total_cons -= records[a][seq[a]]
	while a <= l:
		ss_cons += records[a][seq[a]]
		total_cons += records[a][seq[a]]
		a += 1
	ss_cons /= 7
	while a < len(seq):
		if a in records:
			ds_cons += records[a][seq[a]]
			total_cons += records[a][seq[a]]
		a += 1
	ds_cons /= (a - 21)
	total_cons /= a
	loc_dict['Total conservation'] = total_cons
	loc_dict['Slippery Sequence conservation'] = ss_cons
	loc_dict['Upstream Sequence conservation'] = up_cons
	loc_dict['Downstream Structure conservation'] = ds_cons
	if sd_cons > 0:
		loc_dict['SD conservation'] = sd_cons
	return(loc_dict)