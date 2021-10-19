#!/usr/bin/env python3

import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import subprocess
import pandas as pd
import numpy as np
from collections import OrderedDict
import itertools


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
	md_df[f'{name}_kmer_pos'] = md_df.header.map(pos_map)
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
	pos_name = name + '_kmer_pos'
	pos_dict = OrderedDict()
	# get rows with cycles that have mapped completely to genome
	cycle_df = md[md.header.str.contains('cycle')]
	cycle_idx = cycle_df.where(cycle_df.seq.apply(lambda x: len(x)) == cycle_df[name + '_kmer_pos'].apply(get_len)).dropna(subset=['header']).index
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

def get_mindist(genome_pos, fpos):
	"""
	Get the minimum distance between the starts and the ends of the kmer's genome pos and feature pos
	"""
	#genome pos is always tuple and so is the fpos
	combo = itertools.product(genome_pos, fpos)
	min_dist = 10e9
	for pair in combo:
		min_dist = min(min_dist, abs(pair[0] - pair[1]))
	return min_dist

def maptofeatures(md, feats, name):
    """
    Maps the kmers to the features in the genome; in case of intergenic
	mutations, the nearest feature is mapped instead. Updates the info in the metadata.
	
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
    col = name + '_kmer_pos'
	# {qseqid: (start, end)}
    pos_dict = get_positions_dict(md, name)
    idx = 0
    matched_dict = {}
    for qseqid, pos in pos_dict.items():
        matched = []
        starting_idx = idx
        genome_pos = feats.loc[idx, ['start', 'end']].values

        # map the closest feature, in case the mutation is intergenic
       
        min_dist = get_mindist(genome_pos, pos)
        closest_match = feats.loc[idx, 'locus_tag']

        # keep cycling through genome til you find overlap
        while not hasoverlap(pos, genome_pos):
            try:
                fdist = get_mindist(genome_pos, pos)
                min_dist = min(fdist, min_dist)
                if min_dist == fdist: # new closest feature
                    closest_match = feats.loc[idx, 'locus_tag']

                idx += 1
                genome_pos = feats.loc[idx, ['start', 'end']].values
            except KeyError: # no match at all
                idx = starting_idx
                break

        starting_idx = idx # since sorted, all future matches should be after this

        #found overlap, match can overlap multiple features so check them all
        while hasoverlap(pos, genome_pos):
            try:
              matched.append(feats.loc[idx, 'locus_tag'])
              idx += 1
              genome_pos = feats.loc[idx, ['start', 'end']].values
            except KeyError:
              idx = starting_idx
              break
        if matched:
            lts = ';'.join(matched)
        else:
            lts = closest_match + '(intergenic)'
            
        matched_dict.update({qseqid: lts})
        # need to reset idx and check again since multiple MGE can map to same feats
       	idx = starting_idx
    combined_dict = combine_feat_dict(matched_dict)
    combined_pos_dict = combine_pos_dict(pos_dict)
    md[name + '_feats'] = md['header'].map(combined_dict)
    md[name + '_mutation_pos'] = md['header'].map(combined_pos_dict)

def combine_pos_dict(pdict):
	"""
	Helper function to combine the positions of multiple mutations within same kmer
	"""
	combined_pos_dict = {}
	for key, vals in pdict.items():
		name = key.split('mutation')[0]
		if name in combined_pos_dict:
			combined_pos_dict[name].append(vals)
		else:
			combined_pos_dict[name] = [vals]
	return combined_pos_dict

def combine_feat_dict(matched_dict):
	"""
	Helper function to combine features from multiple mutations in same cycle
	"""
	combined_dict = {}
	for key, vals in matched_dict.items():
		name = key.split('mutation')[0]
		vals = ';'.join(set(vals)) if type(vals) == list else vals
		if name in combined_dict and vals not in combined_dict[name]:
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

def getseqdict(gb):
	seq_dict = {}
	for refseq in SeqIO.parse(gb, 'genbank'):
		for feats in refseq.features:
			try:
				lt = feats.qualifiers['locus_tag'][0]
			except KeyError:
				continue
			else:
				seq_dict.update({lt: str(feats.extract(refseq).seq)})
	return seq_dict
	
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
				continue
			else:
				feats.append([locus_tag, feat.type, int(feat.location.start),
						  int(feat.location.end), feat.location.strand])
	return pd.DataFrame(feats, columns=['locus_tag', 'ftype', 'start', 'end', 'strand'])

def parse_muts_pos(mut):
    """
    Parse mutation positions and find the start and end.
    """
    if not mut:
        return
    
    if abs(mut[0] - mut[1]) == 1:
        return mut[0]
    return np.arange(min(mut), max(mut))

def parse_muts_strand(mut):
    return 1 if mut[0] < mut[1] else -1

def replace_str_index(text,index=0,replacement=''):
	"""
	Helper function to replace the string at the given index
	"""
	if index >= len(text):
		return ValueError('Mutation position greater than the sequence.')
	return f'{text[:index]}{replacement}{text[index+1:]}'

def get_aamutation(codon, ncodon, pos):
	"""
	Helper function to translate codon to aa
	"""
	aa1, aa2 = Seq(codon).translate(), Seq(ncodon).translate()
	return f'{aa2}{pos//3}{aa1}'

def get_snps(mps, sps, fstrand, mt, seq, ftype, pheno):
	"""
	Returns the identity e.g. G23C of SNP mutation within the feature.
	Parameters
	----------
	mps: int
		mutation position within the feature
	sps: int (1, -1)
		strand of the kmer, 1 for pos and -1 for neg
	fstrand: int (1, -1)
		strand of the feature the mutation is mapped to
	mt: str
		btop representation of mutation e.g. 'GC' 'A-' etc.
	seq: str
		sequence of the feature
	ftype: str
		type of feature e.g. CDS, rRNA etc.
	pheno: str, pheno0 or pheno1
		the phenotype of the reference genome where the feature is pulled from
	
	returns
	-------
	sub: str
		mutation information represented by the snp, mapped to the feature position; if
		feature is CDS, the AA substitution is also added.
	"""
	mps = int(mps)
	if fstrand == sps:
		point_mt = seq[mps]
	else:
		mps -= 1
		point_mt = Seq(seq[mps]).reverse_complement()
	check_pos = 0 if pheno == 'pheno0' else 1
	if mt[check_pos] != point_mt:
		return 'No Match'
    
	if point_mt == '-':
		return 'Del' + str(mps)
	elif mt[0] == '-':
		return 'Ins' + str(mps)
	else:
		sub = f'{mt[0]}{mps}{mt[1]}'
        # find the codon string
		if ftype == 'CDS':
			codon_start_adj = mps % 3
			codon_end_adj = 3 - codon_start_adj
			codon = seq[mps - codon_start_adj: mps + codon_end_adj]
			new_codon = replace_str_index(codon, codon_start_adj, point_mt)
			aamutation = get_aamutation(codon, new_codon, mps)
			sub += f'({aamutation})'
		return sub

def get_indel(mps, sps, mt, seq, pheno):
	"""
	Get mutation information for indel or multiple location mutation

	Parameters
	----------
	 mps: int
         mutation position within the feature
     sps: int (1, -1)
         strand of the kmer, 1 for pos and -1 for neg
     mt: str
         btop representation of mutation e.g. 'GC' 'A-' etc.
     seq: str
         sequence of the feature
     ftype: str
         type of feature e.g. CDS, rRNA etc.
     pheno: str, pheno0 or pheno1
         the phenotype of the reference genome where the feature is pulled from
	"""

	current_seq = mt[::2] # sequence of current phenotype
	sub_seq = mt[1::2] # substituations in the other phenotype
    
    # position on btop format mutation string is pheno0pheno1
    # e.g. 37AGTG23 means 'AT' in pheno0 and 'GG' mutation in pheno1
	if pheno == 'pheno1':
		current_seq, sub_seq = sub_seq, current_seq
    
	if len(set(current_seq)) == 1 and set(current_seq) == set('-'): # deletion
		geno = 'Del' + str(min(mps)) + ':' + str(max(mps))
	elif len(set(sub_seq)) == 1 and set(sub_seq) == set('-'): # insertion
		geno = 'Ins' + str(min(mps)) + ':' + str(max(mps))
	else:
		ans = []
		[ans.append(f'{orig}{pos}{sub}') for orig, pos, sub in zip(current_seq, mps, sub_seq)]
		geno = ''.join(ans)
	return geno
						
def get_feature_mutation(cds_md, feats_df, pheno, name):
	"""
	For all mutations mapped to features, get the mutation information
	
	Parameters
	----------
	cds_md: pandas.DataFrame
		metadata of the kmers to be mapped 
	feats_df: pandas.DataFrame
		dataframe containing features the kmers will be mapped to
	pheno: str, pheno0 or pheno1
		phenotype of the reference genome
	name: str
		name of the reference genome
	"""
	mut_map = {}
	for idx, row in cds_md.iterrows():
		mutations = []
		match = feats_df[feats_df.locus_tag == row[name + '_feats']]
		try:
			fstart, fend, fstrand, ftype = match[['start', 'end', 'strand', 'ftype']].values[0]
		except IndexError: #no match
			continue
		start = fstart if fstrand == 1 else fend
		num = 0
		for mt in re.split(r'(\d+)', row.mutations)[1:-1]:
			try:
				int(mt)
			except ValueError:
				mut = parse_muts_pos(row[f'{name}_mutation_pos'][num]) # can be int or range
				sps = parse_muts_strand(row[f'{name}_mutation_pos'][num])
				if type(mut) == np.ndarray: #indels, multiple bases
					if any([fstart <= m <= fend for m in mut]):
						mutations.append(get_indel(mut, sps, mt, 'seq', 'pheno1'))
					else:
						m1, m2 = min(mut), max(mut)
						mutations.append(f'Indel{m1}:{m2}')
				else: # SNPs
					if fstart <= mut <= fend:
						cds_pos = abs(start - mut)
						mutations.append(get_snps(cds_pos, sps,
                                                  fstrand, mt,
                                                  seq_dict[row[f'{name}_feats']],
                                                  ftype, pheno))
					else:
						m1, m2 = mt
						if pheno == 'pheno0':
							m1, m2 = m2, m1
						mutations.append(f'{m1}{mut}{m2}')
				num += 1
		mut_map.update({idx:';'.join(mutations)})
	return mut_map

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
			feats_df = getfeatures(gen)
			seq_dict = getseqdict(gen)
			# find the closest feature
			maptofeatures(md_df, feats_df, name)
			# find the mutation within the feature, including aa mutation if in CDS
			cds_md = md_df[md_df.header.str.contains('cycle')].dropna(subset=[f'{name}_feats'])
			mut_map = get_feature_mutation(cds_md, feats_df, pheno, name)		
			md_df[f'{name}_mutation_map'] = md_df.index.map(mut_map)

	md_df.to_csv('test_md.csv')
		
			
