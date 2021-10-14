#!/usr/bin/env python3

import os
import re
from Bio import SeqIO
import tempfile
import subprocess
import pandas as pd
import numpy as np
from collections import OrderedDict

def makeblastdb(gen, direct):
	"""
	Make blastn database from the genbank file.
	
	Parameters
	----------
	gen: str, path.PATH
		path to the genbank file
	direct: str, path.PATH
		path to directory where the database will be created

	returns
	-------
	fa_file: str
		path to the newly created blasn database
	"""
	fa_file = os.path.join(direct, name + '.fasta')
	with open(fa_file, 'w') as fout:
		for refseq in SeqIO.parse(gen, 'genbank'):
			fout.write(f'>{refseq.id}\n{refseq.seq}\n')
	#dst = os.path.join(direct, os.path.basename(fa_file))
	#copyfile(fa_file, dst)
	db_call = ['makeblastdb', '-dbtype', 'nucl', '-in', fa_file]
	with open(os.devnull, 'w') as null_buffer:
		subprocess.call(db_call, stdout=null_buffer)
	return fa_file 


def runblast(fa, db, name, direct):
	"""
	Run blastn between the reference database and the component fasta file
	
	Parameters
	----------
	fa: str, path.PATH
		path to the component fasta file
	db: str, path.PATH
		path to the database file
	"""
	global outfmt
	outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
	blast_call = ['blastn', '-db', db, '-query', fa, '-max_target_seqs', '1', '-outfmt', outfmt]
	outfile = os.path.join(direct, name + '.txt')
	with open(outfile, 'w') as tout:
		tout.write(subprocess.check_output(blast_call).decode('utf-8'))
	return outfile

def process_blastout(blast_out, pheno):
	"""
	Parse blast results to find the DBGWAS sequences that match the reference genome
	
	Parameters
	----------
	blast_out: str, path.PATH
		path to the blast result (outfmt 6)

	returns
	-------
	blast_df: pandas.DataFrame
		pandas dataframe containing the sequence info for those that mapped to the reference
	"""

	blast_df = pd.read_csv(blast_out, sep='\t', header=None, names=outfmt.split()[1:]) #outfmt from runblast
	blast_df = blast_df[blast_df.qseqid.str.contains(pheno)]
	# sort by bitscore so when duplicated is dropped, the one with highest is kept
	blast_df.sort_values('bitscore', ascending=False, inplace=True)
	blast_df.drop_duplicates(subset=['qseqid'], inplace=True)
	blast_df.reset_index(inplace=True, drop=True)
	keep_idx = np.concatenate((_mge_match(blast_df)[0], _cycles_match(blast_df)[0]))
	blast_df = blast_df.loc[keep_idx]
	return blast_df

def _mge_match(df):
	"returns index of sequences that match MGE criteria"
	return np.where((df.qseqid.str.contains('MGE') &
                    (df.pident >= 80.0) &
                    (df.qcovs >= 80.0)))

def _cycles_match(df):
	"returns index of cycle sequences that match the reference genome"
	return np.where((df.qseqid.str.contains('cycle')) &
            (df.pident == 100.0) &
            (df.qcovs  == 100.0))

def map_blast(blast_df, md_df, name):
	"""
	Add the sequence info from the blast data to the metadata
	
	Parameters
	----------
	blast_df: pandas.DataFrame
		pandas dataframe containing aln info from mapped DBGWAS sequences (outfmt 6)
	md_df: pandas.DataFrame
        pandas dataframe containing the DBGWAS sequence metadata
	"""
	# for the all sequences there, map (start, end) -> 'pos' 
	pos_map = dict(zip(blast_df.set_index('qseqid').index, zip(blast_df['sstart'], blast_df['send'])))
	md_df[f'{name}_pos'] = md_df.header.map(pos_map)
	# for mutations, convert the positions to

def get_len(prange):
	"""Get length of seq from the range"""
	if prange != prange:
		return np.nan
	n1, n2 = prange
	return abs(int(n1) - int(n2)) + 1

def get_mut_pos(muts, genome_pos):
	"""
	Get the genomic position of SNP or indel mutations.

	Parameters
	----------
	muts: str
		btop formatted string showing position of mutation
	genome_pos: tuple
		start and end of the kmer_pos in the genome
	"""
	mut_positions = []
    #p1, p2 = [int(i) for i in genome_pos[1:-1].split(',')]
	p1, p2 = genome_pos
	fwd = p1 < p2
    # blast range is inclusive while python is exclusive
	if fwd:
		p2 += 1
	else:
		p1 -= 1
    #sum or diff depending on fwd or revcomp
	traverse = lambda x,y: x+y if fwd else x-y
	adj = -1 if fwd else 1
	start = p1
	for vals in re.split(r'(\d+)', muts)[1:-1]:
		try:
			vals = int(vals)
		except ValueError: # not int
			if len(vals) % 2 != 0:
				raise ValueError('The mutation info was not an even number.')
			mut_end = int(traverse(start, (len(vals) / 2)) + adj)
			mut_positions.append((start + adj, mut_end))
			start = traverse(start, len(vals) / 2)
		else:
			start = traverse(start, vals)
	return mut_positions

def get_positions_dict(md, name):
	"""
	Generate the mutation positions of all aligned kmers
	
	Parameters
	----------
	md: pd.DataFrame
		pandas dataframe containing the DBGWAS sequence metadata
	name: str
		name of the reference genome to get the position from
	returns
	-------
	pos_dict: dict, {str:tuple}
		dict containing the name of the mutation and its position in genome
	"""
	pos_name = name + '_pos'
	pos_dict = OrderedDict()
	# get rows with cycles that have mapped completely to genome
	cycle_df = md[md.header.str.contains('cycle')]
	cycle_idx = cycle_df.where(cycle_df.seq.apply(lambda x: len(x)) == cycle_df[name + '_pos'].apply(get_len)).dropna(subset=['header']).index
	mge_idx = md[md.header.str.contains('MGE')].index
	pos_df = md.loc[set(cycle_idx).union(mge_idx)]
	pos_df.dropna(subset=[pos_name], inplace=True)
	pos_df.sort_values(pos_name, inplace=True)
	for idx, row in pos_df.iterrows():
		if 'MGE' in row.header:
			pos_dict.update({row.header: row[pos_name]})
		else:	
			mps = get_mut_pos(row.mutations, row[pos_name])
			for mp, ps in enumerate(mps):
				pos_dict.update({row.header + 'mutation' + str(mp): ps})	

	return pos_dict


def mapcyclefeatures(md, feats, name):
	"""
	Maps mutations to features in the genome, updates the info in the metadata
    md: pd.DataFrame
        pandas dataframe containing the DBGWAS sequence metadata
    feats: pd.DataFrame
        pandas dataframe contaning reference genome feature info
    name: string
        name of the reference
    """

	# filter cycles where the length of the seq path doesn't match the length of the region mapped to the genome
	cycle_df = md.where(md.seq.apply(lambda x: len(x)) == md[name + '_pos'].apply(get_len)).dropna(subset=['header'])
	cycle_df = cycle_df[cycle_df.header.str.contains('cycle')]
	
	# for each cycle, find the location(s) of the mutations
    # if the mutations are in the genes then
	
def maptofeatures(md, feats, name):
    """
    Maps MGEs to the features in the genome, updates the info in the metadata
	
	Parameters
	----------
	md: pd.DataFrame
		pandas dataframe containing the DBGWAS sequence metadata
	feats: pd.DataFrame
		pandas dataframe contaning reference genome feature info
	name: string
		name of the reference
    """
	# dict needs to preserve the sorted order of positions
    col = name + '_pos'
	# {qseqid: (start, end)}
    pos_dict = get_positions_dict(md, name)
    idx = 0
    matched_dict = {}
    for qseqid, pos in pos_dict.items():
        matched = []
        starting_idx = idx
        genome_pos = feats.loc[idx, ['start', 'end']].values
        # keep cycling through genome til you find overlap
        while not hasoverlap(pos, genome_pos):
            try:
                idx += 1
                genome_pos = feats.loc[idx, ['start', 'end']].values
            except KeyError: # no match at all
                idx = starting_idx
                break

        starting_idx = idx # since sorted, all future matches should be after this

        #found overlap, match can overlap multiple features
        while hasoverlap(pos, genome_pos):
            try:
              matched.append(feats.loc[idx, 'locus_tag'])
              idx += 1
              genome_pos = feats.loc[idx, ['start', 'end']].values
            except KeyError:
              idx = starting_idx
              break
        if matched:
            matched_dict.update({qseqid: ';'.join(matched)})
        # need to reset idx and check again since multiple MGE can map to same feats
       	idx = starting_idx
    combined_dict = combine_dict(matched_dict)
    md[name + '_feats'] = md['header'].map(combined_dict)


def combine_dict(matched_dict):
	"""
	Helper function to combine features from multiple mutations in same cycle
	"""
	combined_dict = {}
	for key, vals in matched_dict.items():
		name = key.split('mutation')[0]
		vals = ';'.join(vals) if type(vals) == list else vals
		if name in combined_dict:
			combined_dict[name] = combined_dict[name] + ';' + vals
		else:
			combined_dict[name] = vals
	return combined_dict
	
def hasoverlap(kmer_range, gfeat):
	"""
	Checks if there is overlap between two ranges
	
	Parameters
	----------
	kmer_range: tuple
		start and end positions of the compacted kmer
	gfeat: tuple
		start and end positions of the genomic feature
	
	returns
	-------
	overlap: boolean
		whether the two ranges overlap
	"""

	# snp
	if abs(kmer_range[0] - kmer_range[1]) == 1:
		return gfeat[0] <= kmer_range[0] <= gfeat[1]

	a, b = sorted(kmer_range)
	c, d = sorted(gfeat)
		
	return max(a,c) <= min(b,d)
	
	
def getfeatures(gb):
	""" 
	Get all the information about different features in the genbank files
	
	Parameters
	----------
	gb: str, path.PATH
		path to the genbank file
	
	returns
	-------
	feats: pd.DataFrame
		dataframe containing feature information- locus tag, feature type, start, end, and strand
	"""
	feats = []
	for rs in SeqIO.parse(gb, 'genbank'):
		for feat in rs.features:
			if feat.type in ['gene', 'source']:
				continue
			# locus_tag, ftype, start, end
			try:
				locus_tag = feat.qualifiers['locus_tag'][0]
			except KeyError:
				locus_tag = 'Unknown feature'
			feats.append([locus_tag, feat.type, int(feat.location.start),
						  int(feat.location.end), feat.location.strand])
	return pd.DataFrame(feats, columns=['locus_tag', 'ftype', 'start', 'end', 'strand'])
						

	
if __name__ == '__main__':
	
	import argparse
	p = argparse.ArgumentParser(description='Annotate the mutations using reference genomes. The annotations are updated in the input\
								metadata file.')
	p.add_argument('-f','--fasta', help='path to fasta file containing the path and MGE sequences. Generated by component_analysis.py',
					type=str)
	p.add_argument('-m', '--md', help='path to metadata file containing sequnece metadata. Generated by component_analysis.py',
					type=str)
	p.add_argument('-r', '--ref', help='Reference genbank files used to annotate the sequences. Must be reference sequence for pheno0 followed by reference for pheno1',
				    required=True, nargs='+')
	p.add_argument('-p', '--pheno', help='Whether the references passed in -r are from pheno0 or pheno1. Must pass same number of values as number of references e.g. "-p pheno0 pheno1"',
				    required=True, nargs='+') 
	params = vars(p.parse_args())
	
	#if len(params['ref']) != 2:
		#raise ValueError(f'Must provide 2 reference genomes, one for pheno0 and one for pheno1 in that order. {len(params["ref"])} genomes provided instead')
	
	md_df = pd.read_csv(params['md'], index_col=0)
	if 'mutations' not in md_df.columns:
		raise AttributeError('"Mutations" column not found in metadata. Run "sequence_analysis.py" first')
	
	if len(params['pheno']) != len(params['ref']):
		raise ValueError('Different number of reference and phenotypes passed. Recheck number of params passed to "-r" and "-p."')
	with tempfile.TemporaryDirectory() as tempdir:
		for pheno, gen in zip(params['pheno'], params['ref']):
			name = os.path.basename(gen).split('.')[0]
			db_file = makeblastdb(gen, tempdir) #make db of reference genome
			blast_outfile = runblast(params['fasta'], db_file, name, tempdir)
			matched_seq = process_blastout(blast_outfile, pheno)
			map_blast(matched_seq, md_df, name)
			feats = getfeatures(gen)
			maptofeatures(md_df, feats, name)
	print('SAVING')
	md_df.to_csv('test_md.csv')
		
			
