import argparse
import os
import datetime

def split_fasta(input_fasta, output_folder):
	# Splits fasta with genes into folder for each gene
	# Input: fasta with gene
	# Output: directory with folders for each gene sequence
	
	with open(input_fasta) as f:
		combined_fasta = f.read()
		combined_fasta = combined_fasta[1:]
	f.close()
	#initialize file_name variable
	file_name = "init"
	#initialize dictionary with name/ content of each new fasta
	name_content_dict = {}
	genes_list = combined_fasta.split('\n>')
	#(reference[i:i+3] not in stop_codons) and
	gene_lengths = {}
	for gene in genes_list:
		# extract name of gene from header
		#file_name = gene.split(' ')[0]
		header = gene.split('\n')[0]
		seq = ''.join(gene.split('\n')[1:])
		file_name = header.split(':')[3]
		gene_lengths[file_name] = len(seq)
		if len(file_name) == 0:
			continue
		# put gene name and sequence (header included) into dictionary
		name_content_dict[file_name] = '>' + gene
	#create new files/ directories
	os.mkdir(output_folder)
	for name, content in name_content_dict.items():
		os.mkdir('{}/{}'.format(output_folder, name))
		with open('{}/{}/{}.fa'.format(output_folder, name, name), 'w+') as f:
			f.write(content)
			f.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', default=None, type=str, help='Path to fasta')
	parser.add_argument('-O', default=None, type=str, help='output folder')
	args = parser.parse_args()
	split_fasta(args.F, args.O)