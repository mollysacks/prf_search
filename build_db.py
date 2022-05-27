import requests
import csv
import io
import gzip
import argparse
import pandas as pd

def build_db(inp, out):
	df = pd.read_csv(inp, sep='\t')
	total_genomes = len(df)
	i = 0
	genomes = []
	for url in df.ftp_path:
		print(i/total_genomes)
		zipf = url + '/' + url.split('/')[-1] + '_genomic.fna.gz'
		response = requests.get(zipf, stream=True)
		content = str(gzip.decompress((io.BytesIO(response.content).read())))
		content = content.replace(r'\n', '\n')
		genomes.append(content[2:])
		i += 1
		with open(out, 'a') as f:
			f.write('\n'.join(genomes))
		f.close()
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-I', default=None, type=str, help='Path to assembly summary')
	parser.add_argument('-O', default=None, type=str, help='Path to output file')

	args = parser.parse_args()
	build_db(args.I, args.O)

	