#!/usr/bin/env python3

from Bio import Align, Seq, SeqIO
import pandas as pd
import numpy as np
import os
import itertools

# TODO: remove this, it is now depracated in favor of global aln
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

def align_it(seq1, seq2, aligner):
	"""
	Align seq1 and seq2 against one another using global alignment
	
	Parameters
	----------
	seq1/seq2: str
		sequences to be aligned, assumes DNA sequences containing only 'ATCG'
	aligner: SeqIO.Align
		aligner object initialized with the desired penalty values etc.
	
	returns
	-------
	btop: str
		btop formatted string describing the alignment between seq1 and seq2
	"""
	alignments = aligner.align(seq1, seq2)
	if not alignments:
		return
	max_score = np.argmax([aln.score for aln in alignments])
	return get_btop(alignments[max_score])

def get_btop(align):
	"""
	Converts alignment sequence to a more compact btop format.
	
	Parameters
	----------
	aln: str
		alignment string, default output of biopython alignment programs
	
	returns
	-------
	btop: str
		btop formatted string describing the aln string
	"""

	idx = 0
	btop = ''
	seq1, aln, seq2 = str(align).split()
	for name, group in itertools.groupby(aln):
		grouplen = sum(1 for match in group)
		if name == '|':
			idx += grouplen
			btop += str(grouplen)
			continue
		for g in np.arange(grouplen):
			btop += str(seq1[idx]) + str(seq2[idx])
			idx += 1
	return btop

def yield_paths(fa):
	"""
	Yields sequence of paths from same cycle. Different sequence paths in same cycle represent
	different DBGWAS enriched mutations.
	
	Parameters
	----------
	fa: str
		path to fasta file containing path sequences, assumes that different path sequnces
		in same cycle are in adjacent positions in the fasta file
	
	Yields
	------
	seqs: list
		list containing the name of current cycle and sequnce of the paths in that cycle
	Returns
	-------
		returns 0 when EOF is reached
	"""
	parse_fa = SeqIO.parse(fa, 'fasta')
	rs = next(parse_fa)
	cycle = ''
	while(rs):
		if 'MGE' in rs.id or rs.id.split('pheno')[0] == cycle:
			try:
				rs = next(parse_fa)
				continue
			except StopIteration:
				return 0

		# pheno0 of cycle
		id1, seq1 = rs.id, rs.seq
		cycle = rs.id.split('pheno')[0]
		rs = next(parse_fa)
		id2, seq2 = rs.id, rs.seq
		if id2.split('pheno')[0] != cycle:
			print('IDS', id1, id2)
			raise ValueError('Sequence id mismatch; make sure fasta is sorted i.e. >pheno0 followed by >pheno1 or vice-versa')
		# switch so that pheno0 is always id1, seq1
		if 'pheno1' in id1:
			id1, id2 = id2, id1
			seq1, seq2 = seq2, seq1
		yield cycle, seq1, seq2			
		
		if rs is None: #EOF
			return 0

#TODO: remove, depracated
def groupseq(parse_fa, rs, cycle_name, direct):
	"""
	 Helper function, to group sequences that are from component and cycle
	"""
	seqids = []
	seqs = []
	while rs.id.split('pheno')[0] == cycle_name:
		seqids.append(rs.id)
		seqs.append(rs.seq)
		try:
			rs = next(parse_fa)
		except StopIteration:
			rs = None
			break
	return seqids, seqs, rs

def single_cycle_analysis(seqs, aligner):
	"""
	Compares the sequences from the paths of the single cycle and returns something
	
	Parameters
	----------
	seqs: list
		list of tuple, each tuple containing the sequence id and sequence to be compared
	returns
	-------
	mutations:	str
		btop formatted str containing the mutation differences between the two path seqs e.g. '23GT'
		pheno0 is always the reference so in this example pheno0 has 'G' and 
		pheno1 has 'T' in position 24 after 23 matches.
	"""
	cycle, seq1, seq2 = seqs
	return cycle, align_it(seq1, seq2, aligner)

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
	md['mutations'] = md.networkFeature.map(dict(mutations))
	md.to_csv(md_file)

def initialize_aligner():
	"""
	Helper function to initialize aligner and keep main clean
	"""
	aligner = Align.PairwiseAligner()
	aligner.mode = 'global'
	aligner.match_score = 2
	aligner.mismatch_score = -1
	aligner.open_gap_score = -10
	aligner.extend_gap_score = -0.5
	aligner.end_extend_gap_score = -0.5
	aligner.end_open_gap_score = -10
	return aligner
	
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
	aligner = initialize_aligner()
	mutations = [single_cycle_analysis(seq_info, aligner) for seq_info in yield_paths(params['fa_file'])]
	
	update_metadata(mutations, params['md_file'])
	
