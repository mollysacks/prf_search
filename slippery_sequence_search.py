import argparse
import re
import os
import shutil
import json

def search(reference, length, min_struct_len, upstream_len, space, gene_length):
	# searches through sequence for slippery sequences
	# outputs list of tuples with (loca, sequence + L bases)
	segments = list() #initialize list of segments
	i = 26 # start at last base of ninth codon
	while (i <= gene_length): #len(reference) - 7): # go until last 7 bases
		if re.match(r'([ATCG])\1{2}([AT])\2{2}[ATCG]', reference[i:i+7]): # check if slippery sequence
			if len(reference) - i > length: # if there are > length bases left
				new_seq = reference[i - (upstream_len + 12): i + space + length] # include an extra 12 bases in case hmmer cuts some off
			segments.append((i, new_seq))
		i = i + 3 # go to next codon
	return segments

def realign_newlines(seq):
	# write new sequence with properly spaced \n
	new_seq = ''
	i = 0
	while i < len(seq):
		j = i + 79
		if j > len(seq):
			new_seq = new_seq + seq[i:]
		else:
			new_seq = new_seq + seq[i:j] + '\n'
		i += 79
	return new_seq

def find_slippery_sequences(rootdir, length, min_struct_len, upstream_len, space):
	# rootdir contains a directory for each cDNA sequence in input file
	# Each of these folders contains a fasta w/ just that cDNA sequence
	with open('gene_lengths_e_coli') as f:
		glength = f.read()
		glength = glength.replace("'", '"')
		glength = json.loads(glength)
		f.close()
	for directory in os.listdir(rootdir):
		for fasta in os.listdir(os.path.join(rootdir, directory)):
			with open(os.path.join(rootdir, directory, fasta)) as f:
				reference = f.read()
			f.close()
			lines = reference.split("\n")
			header = lines[0]
			gene = header.split(':')[3]
			sequence = ''.join(lines[1:])
			if not re.match(r'.TG', sequence[0:3]):
				continue
			if gene not in glength:
				glength[gene] = len(sequence)
			segments = search(sequence, length, min_struct_len, upstream_len, space, glength[gene])
			if segments == []:
				shutil.rmtree(os.path.join(rootdir, directory))
			for loc, seq in segments:
				# realign newlines to improve readability
				seq = realign_newlines(seq)
				os.mkdir(os.path.join(rootdir, directory, '{}-{}'.format(directory, loc)))
				with open(os.path.join(rootdir, directory, '{}-{}'.format(directory, loc),'{}-{}.fa'.format(directory, loc)), "w+") as f:
					f.write(header + "\n")
					f.write(seq)
					f.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-D', default=None, type=str, help='directory w cDNA folders')
	parser.add_argument('-L', default=None, type=int, help='length of sequence after slippery to include in output')
	parser.add_argument('-M', default=18, type=int, help='min length of downstream sequence to include in output')
	parser.add_argument('-U', default=90, type=int, help='max length of upstream segment to search for SD and nascent polypeptide')
	parser.add_argument('-B', default=12, type=int, help='distance between slippery start and potential downstream structure')
	args = parser.parse_args()
	find_slippery_sequences(args.D, args.L, args.M, args.U, args.B)


