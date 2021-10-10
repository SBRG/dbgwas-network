#!/usr/bin/env python3

import networkx as nx
from networkx.exception import NetworkXNoPath
from networkx.algorithms.shortest_paths.generic import shortest_path
from networkx.readwrite import json_graph

import json
import numpy as np
from Bio.Seq import Seq
import itertools
import pandas as pd
import glob
import os 
np.seterr('raise')

def cytoscape_graph(data, mincount, attrs=None, name="name", ident="id"):
    """
    NOTE: Adapted from NetworkX, modified to deal with discrepency in json created by
    DBGWAS method. 
    Create a NetworkX graph from a dictionary in cytoscape JSON format.

    Parameters
    ----------
    data : dict
        A dictionary of data conforming to cytoscape JSON format.
    mincount: int
        Int describing minimum number of genomes the node must appear in. Nodes
        with count lower than this will be excluded from the graph.
    attrs : dict or None (default=None)
        A dictionary containing the keys 'name' and 'ident' which are mapped to
        the 'name' and 'id' node elements in cyjs format. All other keys are
        ignored. Default is `None` which results in the default mapping
        ``dict(name="name", ident="id")``.

        .. deprecated:: 2.6

           The `attrs` keyword argument will be replaced with `name` and
           `ident` in networkx 3.0

    name : string
        A string which is mapped to the 'name' node element in cyjs format.
        Must not have the same value as `ident`.
    ident : string
        A string which is mapped to the 'id' node element in cyjs format.
        Must not have the same value as `name`.

    Returns
    -------
    graph : a NetworkX graph instance
        The `graph` can be an instance of `Graph`, `DiGraph`, `MultiGraph`, or
        `MultiDiGraph` depending on the input data.

    Raises
    ------
    NetworkXError
        If the `name` and `ident` attributes are identical.

    See Also
    --------
    cytoscape_data: convert a NetworkX graph to a dict in cyjs format

    References
    ----------
    .. [1] Cytoscape user's manual:
       http://manual.cytoscape.org/en/stable/index.html

    Examples
    --------
    >>> data_dict = {
    ...     'data': [],
    ...     'directed': False,
    ...     'multigraph': False,
    ...     'elements': {'nodes': [{'data': {'id': '0', 'value': 0, 'name': '0'}},
    ...       {'data': {'id': '1', 'value': 1, 'name': '1'}}],
    ...      'edges': [{'data': {'source': 0, 'target': 1}}]}
    ... }
    >>> G = nx.cytoscape_graph(data_dict)
    >>> G.name
    ''
    >>> G.nodes()
    NodeView((0, 1))
    >>> G.nodes(data=True)[0]
    {'id': '0', 'value': 0, 'name': '0'}
    >>> G.edges(data=True)
    EdgeDataView([(0, 1, {'source': 0, 'target': 1})])
    """
    # ------ TODO: Remove between the lines in 3.0 ----- #
    if attrs is not None:
        import warnings

        msg = (
            "\nThe `attrs` keyword argument of cytoscape_data is deprecated\n"
            "and will be removed in networkx 3.0.\n"
            "It is replaced with explicit `name` and `ident` keyword\n"
            "arguments.\n"
            "To make this warning go away and ensure usage is forward\n"
            "compatible, replace `attrs` with `name` and `ident`,\n"
            "for example:\n\n"
            "   >>> cytoscape_data(G, attrs={'name': 'foo', 'ident': 'bar'})\n\n"
            "should instead be written as\n\n"
            "   >>> cytoscape_data(G, name='foo', ident='bar')\n\n"
            "The default values of 'name' and 'id' will not change."
        )
        warnings.warn(msg, DeprecationWarning, stacklevel=2)

        name = attrs["name"]
        ident = attrs["ident"]
    # -------------------------------------------------- #

    if name == ident:
        raise nx.NetworkXError("name and ident must be different.")

    multigraph = data.get("multigraph")
#     directed = data.get("directed")
    if multigraph:
        graph = nx.MultiGraph()
    else:
        graph = nx.Graph()
    graph.graph = dict(data)
#     graph = graph.to_directed()
    excluded_nodes = []
    for d in data["elements"]["nodes"]:
        # do not add rare nodes
        if int(d['data']['total']) < mincount:
            excluded_nodes.append(d['data'].get(ident))
            continue
        node_data = d["data"].copy()
        node = d["data"]
        if d["data"].get(name):
            node_data[name] = d["data"].get(name)
        if d["data"].get(ident):
            node_data[ident] = d["data"].get(ident)

        graph.add_node(node['id'])
        graph.nodes[node['id']].update(node_data)

    for d in data["elements"]["edges"]:
        edge_data = d["data"].copy()
        sour = d["data"]["source"]
        targ = d["data"]["target"]
        # dont draw edge if one of the node was excluded
        if sour in excluded_nodes or targ in excluded_nodes:
            continue
        if multigraph:
            key = d["data"].get("key", 0)
            graph.add_edge(sour, targ, key=key)
            graph.edges[sour, targ, key].update(edge_data)
        else:
            graph.add_edge(sour, targ)
            graph.edges[sour, targ].update(edge_data)
    return graph


def read_json_file(filename, **kwargs):
    with open(filename) as f:
        js_graph = json.load(f)
    return cytoscape_graph(js_graph, **kwargs)

def get_endnodes(csgraph, graph):
    """Find the 'endnodes' of the cycle i.e. nodes that connect the cycle to the larger graph.
    Assumes that the endnodes of a cycle have the largest total value.
    Parameters
    ----------
    csgraph: networkx.Graph
        the subgraph containing only the cycle    
    Returns
    -------
    nodes: tuple
        size 2 tuple of endnode ids.
    """
    # nodeends are where bifurcations from mutation start/end so we assume that node
    # will appear in most genome compared to other nodes in the cycle and they will have degree > 2
    high_deg_nodes = [i for i in csgraph.nodes if graph.degree[i] > 2]
    node_totals = np.array([int(csgraph.nodes[i]['total']) for i in high_deg_nodes])
    end1, end2 = np.argpartition(node_totals, -2)[-2:]
    return high_deg_nodes[end1], high_deg_nodes[end2]


def get_paths(csgraph, endnodes):
    """
    Returns the two paths through the cycle starting and ending at the same endnodes.
    csgraph: networkx.Graph
        the subgraph containing only the cycle
    endnodes: networkx.nodes
        nodes that signify the 'ends' of the cycle i.e. nodes that connect the 
        cycle to the larger graph.
    Returns
    -------
    nodes: tuple
        size 2 tuple of list containing nodes that form the two paths around the cycle,
        starting and ending at the same endnodes.
    """
    links = list(nx.traversal.dfs_preorder_nodes(csgraph, endnodes[0]))
    endidx = links.index(endnodes[1])
    p1, p2 = links[:endidx + 1], links[endidx:]
    p2.append(endnodes[0])
    
    if p1[0] != p2[-1] or p1[-1] != p2[0]:
        raise ValueError('The ends of the two paths do not match! Something went wrong.')
    return p1, p2[::-1]

def get_nodeconnections(csgraph, path):
    """
    Returns the sequence made from the nodes in the path
    """
    seq = [Seq(csgraph.nodes[p]['name']) for p in path]
    revcomp = [s.reverse_complement() for s in seq]
    all_nodes = [(nodes, {'sequence': sq, 'orientation': 'fwd'}) 
                 for nodes, sq in zip(path, seq)]
    all_nodes.extend([(nodes + '_rc', {'sequence': sq, 'orientation': 'revcomp'})
                      for nodes, sq in zip(path, revcomp)])
    vitree = nx.Graph()
    vitree = vitree.to_directed()
    vitree.add_nodes_from(all_nodes)
    # comapare each node with the subsequent ones to get match
    _ = [getmatch(e1, e2, vitree) for e1, e2 in list(zip(path, path[1:]))]
    
    return vitree

# def endOverlap(a, b):
#     for i in range(0, len(a)):
#         if b.startswith(a[-i:]):
#             return i
#     return 0

def endOverlap(a, b):
    for i in range(0, len(a)):
        if b.startswith(a[i:]):
            return len(a) - i
    return 0

def getmatch(node1, node2, G):
    "Compares whether the ends of the two nodes have overlap"
    # check the k-1 mer for overlap
    nproduct = itertools.product([node1, node1 + '_rc'], [node2, node2 + '_rc'])
    edges = []
    for n1, n2 in nproduct:
        overlap = endOverlap(G.nodes[n1]['sequence'],
                             G.nodes[n2]['sequence'])
        # TODO: need to check if there is a match with the reverse_comp,
        # if so, find the one with the larger overlap
        if not overlap:
            continue
        try:
            # edge with revcomp exists and its overlap is smaller
            check2 = n2[:-3] if '_rc' in n2 else n2 + '_rc'
            if G.get_edge_data(n1, check2)['overlap'] < overlap:
                edges.append((n1, n2, {'overlap': overlap}))
        except TypeError as e: #edge doesn't already exist
            edges.append((n1, n2, {'overlap': overlap}))
    G.add_edges_from(edges)

def concat_paths(path, vitree, name, endorient=None):
    # find all paths through the start and end
    pstarts = [path[0], path[0] + '_rc']
    pends = [path[-1], path[-1] + '_rc']
    product = itertools.product(pstarts, pends)
    spaths = [nx.simple_paths.all_simple_paths(vitree, source, target, cutoff=len(path))
             for source, target in product]
    spaths = [list(i) for i in spaths]
    spaths = [i[0] for i in spaths if i] # get rid of some nested lists
    
    if not any(spaths):
        print(f'\tNo paths found for {name}; check manually')
        return []
    # for each path find the one with the highest overlap
    if len(spaths) == 1:
        return spaths[0]
    # more than one possible path
	# choose the where the orientation of endnodes matches that of the already calculated path through the cycle
    if endorient:
        chosen_path = [i for i in spaths if i[0] == endorient[0] and i[-1] == endorient[-1]]
        if chosen_path:
            #print(chosen_path[0])
            return chosen_path[0]
	# choose path with the highest overlap between sequences
    overlaps = [get_total_overlap(vitree, spath) for spath in spaths]
    return spaths[overlaps.index(max(overlaps))]

def get_total_overlap(vitree, spath):
    return sum([vitree.edges[(n1,n2)]['overlap'] for n1, n2 in zip(spath, spath[1:])])

def get_seqfrompaths(vitree, vitree_path):
    if len(vitree_path) == 0:
        return
    seq = vitree.nodes[vitree_path[0]]['sequence']
    for n1, n2 in zip(vitree_path, vitree_path[1:]):
        seq += concatseq(vitree, n1, n2)
    return seq

def concatseq(vitree, n1, n2):
    """
    Concats sequences from the nodes n1 and n2, removing any overlap from n2. 
    ----------
    vitree: networkx.Graph
        graph that contains the target nodes n1, n2
    n1, n2: networkx.Graph.node
        nodes to concatenate the sequence from
    returns
    -------
    cseq: str
        returns sequnce from n2 with overlapping sequence from n1 removed 
    """
    
    overlap = vitree.get_edge_data(n1, n2)['overlap']
    return vitree.nodes[n2]['sequence'][overlap:]

def hassignodes(graph, path):
    """
    Returns whether there are any 'significant' node i.e. node enriched in DBGWAS
    within the path
    ----------
    graph: networkx.Graph
        graph that contains the nodes in the path
    path: list
        list of nodes in the path
    returns
    -------
    sig: boolean
        returns True if any of the nodes in the path are 'significant'
    """
    for nd in path:
        try:
            if float(graph.nodes[nd]['qValue']) < 0.05:
                return True
        except ValueError as e:
            continue
    #return False

def pathtopheno(graph, path):
    """
    Returns the more frequent phenotype in the path. Assumes binary phenotypes
    labeled 'pheno0' or 'pheno1'
    ----------
    graph: networkx.Graph
        graph that contains the target node
    path: list
        list of nodes in the path
    returns
    -------
    pheno: str
        name of the more frequent phenotype; either 'pheno0' or 'pheno1'
    """
    p0_freq = [calc_node_frequency(graph, i, 'pheno0') for i in path]
    p1_freq = [calc_node_frequency(graph, i, 'pheno1') for i in path] 
    p0mean, p1mean = np.mean(p0_freq), np.mean(p1_freq)
    return 'pheno0' if max(p0mean, p1mean) == p0mean else 'pheno1'


def calc_node_frequency(graph, node, pheno):
    """
    Helper function. Calculate frequency of the node in the given phenotype
    Parameters
    ----------
    graph: networkx.Graph
        graph that contains the target node
    node: networkx.Graph.node
        node to calculate pheno frequency for
    pheno: str
        the phenotype to calculate frequency for; either pheno0 or pheno1 in DBGWAS
    
    returns
    -------
    frequency: float
        frequency of the phenotype
    """
    num, denom = graph.nodes[node][pheno].split('/')
    return float(num) / float(denom)

def write_to_fasta(header, seq, fout):
    """
    Writes the fasta sequnces to file
    Parameters
    ----------
    header: str
        header for the fasta sequence (w/o the starting '>')
    seq: str
        DNA sequence for the feature
    fout: file handle
        handle to write to file
    """
    header = '>' + header
    fout.write(f'{header}\n{seq}\n')
    
def write_metadata(metadata, md_out):
    """
    Writes the metadata to the file
    Parameters
    ----------
    metadata: list
        list of lists with each inner list containing the metadata for each sequnce
        in the fasta file. 
    md_out: str, os.PATH
        path to the file where the metadata will be saved
    """
    cnames = ['header', 'seq', 'nodes']
    pd.DataFrame(metadata, columns=cnames).to_csv(md_out)

def check_mge(comp_graph):
    """
    Checks whether the graph have mobile genetic element (MGE). MGE are defined
    as components where all significant nodes have higher frequency in the same 
    class e.g. all sig nodes are found in pheno0 only. MGE as detected here are often 
    (but not always) part of larger elements and are seldom whole MGE themselves. 
    Parameters
    ----------
    comp_graph: networkx.Graph
        component graph generated by DBGWAS to check for MGE
        
    returns
    -------
    mge: boolean/str
        if MGE is detected, returns the associated phenotype (pheno0 or pheno1) else 
        returns False
    """
    signodes = [nd for nd in comp_graph.nodes if hassignodes(comp_graph, [nd])]
    #if len(signodes) < 2:
        #return False
    p0_freq = [calc_node_frequency(comp_graph, nd, 'pheno0') for nd in signodes]
    p1_freq = [calc_node_frequency(comp_graph, nd, 'pheno1') for nd in signodes]
    check_max = []
    for i in np.arange(len(signodes)):
        check_max.append('pheno0' if p0_freq[i] > p1_freq[i] else 'pheno1')
    if len(set(check_max)) == 1: # same pheno w/ higher frequency in all signodes
        return check_max[0]
    return False

def get_longest_path(comp_graph):
    """
    Returns the longest simple path through the graph that contains atleast one sig node
    Useful for figuring out the nodes to use for MGE sequence.
    Parameters
    ----------
    comp_graph: networkx.Graph
        component graph generated by DBGWAS to check for MGE
        
    returns
    -------
    longest_path: list
        list of nodes in order that form the longest simple path
    """
    node_combinations = itertools.combinations(comp_graph.nodes, 2)
    all_paths_gen = []
    for source, dest in node_combinations:
        try:
            all_paths_gen.append(shortest_path(comp_graph,source, dest))
        except NetworkXNoPath: #nodes not connected once low maf nodes taken out
            continue
    all_paths = [list(i) for i in all_paths_gen]
    signode = False
    while not signode: # longest path must also have at least one sig node
        longest_path = max(all_paths, key=lambda x: len(x))
        signode = hassignodes(comp_graph, longest_path)
        all_paths.remove(longest_path)
    return longest_path

def process_feature(comp_graph, path, pheno, fout, header, endorient=None):
    feat_subgraph = nx.subgraph(comp_graph, path)
    feat_pathgraph = get_nodeconnections(feat_subgraph, path)
    feat_seqpath = concat_paths(path, feat_pathgraph, name=header, endorient=endorient)
    if not feat_seqpath:
        return 0, 0
    feat_seq = get_seqfrompaths(feat_pathgraph, feat_seqpath)
    write_to_fasta(header, feat_seq, fout)
    return feat_seq, feat_seqpath

def jsontoseq(json_dir, tgen, minmaf=0.1, fasta_out='component_seqs.fa',
             md_out='components_md.csv'):
    """
    Create sequence fasta files of all 'genetic events' from all the json network files in
    the directory. The fasta file header can be used to identify the sequence origin and will
    have the format component<X>cycle<Y>pheno<[0,1]>path<Z>. Here a 'cycle' represents a genetic
    event such as SNP or indel that causes bifurcation in the graph while 'path' represents the
    different paths through different bifurcation. Each path is formed by nodes whose sequences
    are associated with a phenotype. The metadata is stored in the pandas file.
    
    Parameter
    ---------
    json_dir: str, os.PATH
        path to the directory that contains the cytoscape json files created by DBGWAS 
    tgen: int
        total number of genomes used in the DBGWAS analysis
    minmaf: float, default 0.1
        minimum mean allele frequency. All nodes with allele frequency lower than this will
        be removed. Can be any float between 0 and 1. Setting it at 0 may result in large
        number of sequences as each individual unique sequence is recorded.
    fasta_out: str, os.PATH, default 'genetic_kmers.fa'
        path to file where the fasta formatted sequences will be written
    md_out: str, os.PATH, default 'genetic_md.csv'
        path to file where the complete metadata for all the sequences will be written
    """
    metadata = []
    with open(fasta_out, 'w') as fout:
        for jsonf in glob.glob(json_dir + '*.json'):
            
            comp = os.path.basename(jsonf).replace('.json', '').split('_')[-1]
            print('Processing component', comp)
            comp_graph = read_json_file(jsonf, mincount=int(minmaf * tgen))
            cycles = nx.algorithms.cycle_basis(comp_graph)
            print(f'\t{len(cycles)} cycle(s) detected')
            
            mge = check_mge(comp_graph)
            if mge: # change to walrus operator on python 3.8
                print(f'\tPossible {mge} associated MGE detected.')
                header = f'component{comp}MGE{mge}'
                mge_path = get_longest_path(comp_graph)
                seq, truepath = process_feature(comp_graph, mge_path, mge, fout, header)
                if seq and truepath:
                    metadata.append([header, seq, ';'.join(truepath)])
                
            for cnum, cyc in enumerate(cycles):
                if not hassignodes(comp_graph, cyc):
                    print(f'\tNo SigNodes; skipping {cyc}')
                    continue
                csgraph = nx.subgraph(comp_graph, cyc) # get the graph of cycles
                endnodes = get_endnodes(csgraph, comp_graph) # get the 'ends' of the cycle	
                # get all the paths across the cycle, one for each mutation
                p1, p2 = get_paths(csgraph, endnodes)
                if len(p1) < 3 or len(p2) < 3: #TODO: fix this, some graphs have mutliple nodes with deg >3,
                    continue
                for cspath in p1, p2:    
                    phenotype = pathtopheno(csgraph, cspath[1:-1])
                    header = f'component{comp}cycle{cnum}{phenotype}'
                    endorient = (truepath[0], truepath[-1]) if cspath == p2 else None
                    seq, truepath = process_feature(comp_graph, cspath, phenotype, fout, header, endorient)
                    if seq and truepath:
                        metadata.append([header, seq, ';'.join(truepath)])
    write_metadata(metadata, md_out)

if __name__ == '__main__':
    import argparse
    
    p = argparse.ArgumentParser(description='Analyze De Bruijn graph of each component to generate fasta sequnce and metadata files')
    
    p.add_argument('json_dir', help='path to directory where the network json files created by \'download_json.py\' is saved.', type=str)
    p.add_argument('tgen', help='total number of genomes used in DBGWAS.', type=int)
    p.add_argument('-m', '--minmaf', help='minimum allele frequency. All nodes with allele frequency lower than this will be discarded', type=float, default=0.1)
    p.add_argument('-f', '--fasta_out', help='path to file where the fasta seq will be saved',
                   type=str, default='component_seqs.fa')
    p.add_argument('-d', '--md_out', help='path to file where the metadata will be saved',
                   type=str, default='components_md.csv')
    params = vars(p.parse_args())
    jsontoseq(**params)
