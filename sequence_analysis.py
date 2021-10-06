#!/usr/bin/env python3

from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import tempfile
import os
import subprocess

def blast_it(query_fa, feat_fasta, dna=True):
	"""
	Make BLASTdb and run BLAST
	
	Parameters
	----------
	query_fa: str, path.PATH
		path to the query fasta file
	feat_fasta: str, os.PATH
		path to the feature fasta file
	blastout: str, path.PATH
		path to file where BLAST output will be saved
	dna: bool, default True
		whether the query and fasta are DNA sequences
	"""
	
	dtype = 'nucl' if dna else 'prot'
	makeblastdb = f"makeblastdb -in {feat_fasta} -dbtype {dtype}"
	with open(os.devnull, 'w') as null_buffer:
		subprocess.call(makeblastdb.split(), stdout=null_buffer)
	blastype = 'blastn' if dna else 'blastp'
	blast_call = [blastype, '-query', query_fa, '-db', feat_fasta, '-max_hsps', '1',  '-outfmt', "6 btop"]
	blast_res = subprocess.check_output(blast_call).decode('utf-8').rstrip()
	if blast_res:
		return blast_res
	# no hits
	blast_call.extend(['-word_size', '5'])
	return subprocess.check_output(blast_call).decode('utf-8').rstrip() 


def yield_paths(fa, direct):
	"""
	Yields sequence of paths from same cycle. Different sequence paths in same cycle represent
	different DBGWAS enriched mutations.
	
	Parameters
	----------
	fa: str
		path to fasta file containing path sequences, assumes that different path sequnces
		in same cycle are in adjacent positions in the fasta file
	direct: str
		path to directory where the output file will be saved
	
	Yields
	------
	gene: str
		name of the current cycle e.g. 'component15cycle3'
	Returns
	-------
		returns 0 when at EOF
	"""
	parse_fa = SeqIO.parse(fa, 'fasta')
	rs = next(parse_fa)
	
	while(rs):
		if 'MGE' in rs.id:
			try:
				rs = next(parse_fa)
				continue
			except StopIteration:
				return 0

		cycle_name = rs.id.split('pheno')[0]
		rs = groupseq(parse_fa, rs, cycle_name, direct)
		yield cycle_name			
		
		if rs is None: #EOF
			return 0

def groupseq(parse_fa, rs, cycle_name, direct):
	fa_locs = []
	while rs.id.split('pheno')[0] == cycle_name:
		fa_loc = os.path.join(direct, rs.id + '.fa')
		with open(fa_loc, 'w') as fout:
			fout.write(f'>{rs.id}\n{rs.seq}\n')
		
		try:
			rs = next(parse_fa)
		except StopIteration:
			rs = None
			break
	return rs

def single_cycle_analysis(cycle, direct, dna=True):
	"""
	Compares the sequences from the paths of the single cycle and returns something
	
	Parameters
	----------
	cycle: str
		name of the current cycle to process
	direct: str, path.PATH
		path where the fasta files are stored
	
	returns
	-------
	mutations:	str
		str containing the mutation differences between the two path seqs e.g. 'G23T.'
		pheno0 is always the reference so in this example pheno0 has 'G' in position 23.
	"""
	# makeblastdb and blast it
	# parse the output and store it in some kind of string 	
	p0_fa = os.path.join(direct, cycle + 'pheno0.fa')
	p1_fa = p0_fa.replace('pheno0.fa', 'pheno1.fa')
	return cycle, blast_it(p0_fa, p1_fa, dna)

def update_metadata(mutations, md_file):
	"""
	Updates the metadata file with the new mutation information
	
	Parameters
	----------
	mutations: list
		list of tuples with each tuple containing the cycle name and the mutations in the paths of the cycle
	md_file: str, path.PATH
		path to the metadata file
	"""
	md = pd.read_csv(md_file, index_col=0)
	md['networkFeature'] = md.header.str.split('pheno', expand=True)[0]
	
	md['mutations'] = md.cycle.map(dict(mutations))
	md.to_csv(md_file)
	
if __name__ == '__main__':
	import argparse
	p = argparse.ArgumentParser('Define mutations for each of the cycle in DBGWAS output')
	p.add_argument('fa_file', help='Path to fasta file containing the sequences from DBGWAS output; created by component_analysis.py',
					default='component_seqs.fa', type=str)
	p.add_argument('md_file', help='Path to metadata file with info on cycles and paths from DBGWAs; created by component_analysis.py',
					default='components_md.csv', type=str)
	p.add_argument('--prots', help='Whether the fasta file contains protein sequence', action='store_false')
	#TODO: Add option to run prots seq, need changes in the single_cycle_analysis function
	params = vars(p.parse_args())
	with tempfile.TemporaryDirectory() as tmpdir:
		mutations = [single_cycle_analysis(cycle, tmpdir, dna=params['prots']) for cycle in yield_paths(params['fa_file'], tmpdir)]

	update_metadata(mutations, params['md_file'])
	
