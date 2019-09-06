#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates criteria of
DiVincenzo and Smolin (cond-mat/9409111). 

The unique gates are iterated by the following steps:
    1) iterate the non-isomorphic graph (this step doesn't implement parallelism)
    2) place edge ordering from 1). At this step, the relabelling qubit has already implemented.
    3) apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping 


MAIN USAGE: 

    from unique2net import unique2net

    unique2net(nqubit, network_length)

"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "3.0.0"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"



#standard libraries
from itertools import product, combinations, permutations 
from time import time
from multiprocessing import Process, Value
from numpy import array_split
from subprocess import run 

import pygraphviz as pgv
import networkx as nx
import json


## bitwise operations CONVENTION: LSB ##
def getbit(num, k):
    """
    get the k-th bit of num, where num has length nqubit

    :num: int, the identified number 
    :k: int, the position

    return 0 or 1
    """
    return (num & (1 << k)) >> k



def swapbits(p1, p2, num):
    """
    Swap p1-th bit with p2-th bit within number num

    :p1:,:p2: int, the swapped position
    :num: int, the number 

    return int"
    Shuffle the binary of num with some permutation 

    :num: int, the number
    :permutaion: tup
    """
    b1 = getbit(num, p1) 
    b2 = getbit(num, p2) 
   
    # XOR the two 
    xor = b1^b2 
   
    # get a mask xor with position 00xor000xor00 
    #                                p1    p2
    xor = (xor << p1) | (xor << p2) 
   
    # or swap with the original number 
    # if xor=0, number is the same, otherwise 0->1 1->0
    res = num ^ xor 
   
    return res 


def get_pos_ones(num):
    """
    Get positions of 1 in num within binary format.
        example: 6 = 110, returns (1,2)
    it is equivalent with extracting the powers of nonzero elements

    :num: int, the number
    """
    max_len = len(bin(num)) - 2   #string format 0bxxxx
    lnum = [int(bool(num & (1<<n))) for n in range(max_len)] #turn into list
    poss = [i for i,e in enumerate(lnum) if e==1] #extract positions

    return tuple(poss)


def shuffle_bits(num, permutation):
    """
    Shuffle the binary of num with some position permutation  

    :num: int, the number
    :permutaion: tuple(int), the final permutation
    """
    num2 = 0  
    for i, new_pos in enumerate(permutation):
        # updating from right to left 0xxx | j000
        num2 =  num2 | (getbit(num, new_pos) << i)

    return num2




## obtaining equivalent networks ##

def equiv_time_reversal(gate_net):
    """you've been a post for a while 
    Equivalent network by reversing the gates order, basically U^{\dagger}

    :gate_net: tuple(int), a configuration of gate network 

    return tuple
    """
    return gate_net[::-1]   #reverse net


def equiv_bit_permutations(gate_net, nqubits):
    """
    Equivalent networks by bit-relabelling, basically performing identical swaps in
    the beginning and the end, equivalent with simple permutation,
    including identity

    :gate_net: tuple(int), a configuration of gate network 
    :nqubits: int, the number of qubits

    return tuple
    """
    return [tuple(shuffle_bits(g, p) for g in gate_net) 
            for p in permutations(range(nqubits))]

        
def equiv_conjugation_by_swapping(gate_net):
    """
    Equivalent network that is equivalent by swapping conjugation.
    Choose the a gates sandwiched by identical gates, then swap it.

    :gate_net: tuple(int), a configuration of gate network 

    return list(tuple)
    """
    result = []
    for j, g in enumerate(gate_net[1:-1]) : #inside big sandwich
        mnet = list(gate_net)  #get a mutable object
        g1, g2 = gate_net[j], gate_net[j+2]  #small sandwich g1|g|g2
        if g1 == g2 :
            poss = get_pos_ones(g1)
            g_swapped = swapbits(*poss, g) 
            mnet[j+1] = g_swapped

            if tuple(mnet)!= gate_net : 
                result.append(tuple(mnet))

    return result


def equiv_DS(nqubit, gate_net, **kwargs):
    """
    Get a set of all equivalent networks based on
    DiVincenzo and Smolin equivalent networks

    :nqubit: int, the number of qubits used
    :gate_list: list(tuple|bool), the list of gate networks

    kwargs arguments : 
        :time_reversal:boolean=False
        :bit_permutations:boolean=True
        :conjugation_by_swap:boolean=True
    """
    opt={'time_reversal':False, 'bit_permutations':True, 'conjugation_by_swap':True} 
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]
    equivs = set()
    if opt['time_reversal']:
        equivs = equivs.union([equiv_time_reversal(gate_net)])
    if opt['bit_permutations']:
        equivs = equivs.union(equiv_conjugation_by_swapping(gate_net))
    if opt['conjugation_by_swap']:
        equivs = equivs.union(equiv_bit_permutations(gate_net, nqubit))
    
    return equivs




## graphs-related methods
def __is_isomorphic_to(G_test, G_list):
    """
    Check if G_test is isomorphic to one of the graph in G_list

    :G_test: nx.graph, the graph to be tested
    :G_list: list(nx.graph), the list of graph to be compared to
    """
    for G in G_list : 
        if nx.is_isomorphic(G,G_test):
            return True
    return False


def __more_three_multiedges(G):
    """
    Check if graph has more than 3 multiple edges

    :G: nx.graph, the tested graph
    """
    L = list(G.edges())
    for edge in L : 
        if L.count(edge) > 3 :
            return True
    
    return False


def list_non_iso_graphs(nqubit, net_depth, **kwargs):
    """
    List the non-isomorphic graphs, where no more that three edges between two nodes.

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks

    kwargs
        :path_json: str, path that store a list of edges of non-isomorphic graphs in json
                    format. If this file presents, the iteration of generating new graphs will be started there.
        :draw_graphs: boolean=True, to draw the produced graphs
        :save_edges: boolean=True, save graph as a collection of edges
        :dirpath: str=out-[n]qubit-depth[d], directory path to store outputs

    """
    #optional arguments
    opt = {'path_json': False, 'draw_graphs': True, 'save_edges': True,i
            'dirpath':'out-'+str(nqubit)+'qubit-depth'+str(net_depth)}
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]

    cphases = list(combinations(range(nqubit),2)) 
    GNIsom = dict() #Non-Isomorphic graphs

    # loading stuff
    if opt['path_json'] : 
        with open(opt['path_json']) as iff :
            list_edges = json.load(iff)
        start_depth = len(list_edges[0])+1
        GNIsom[start_depth-1] = [nx.MultiGraph(edges) for edges in list_edges]
        print("iteration starts from %i edges..."%(start_depth-1))
    else : 
        GNIsom[0]=[nx.MultiGraph()]i
        GNIsom[0][0].add_nodes_from(range(nqubit))
        GNIsom[1]=[nx.MultiGraph()]
        GNIsom[1][0].add_edge(*cphases[0])
        start_depth = 2
        print("iteration starts from 1 edge...")

    # the main part --- iteration of finding non-isomporphic graphs starts here
    for depth in range(start_depth, net_depth+1):
        GNIsom[depth]=list()
        for G_old in GNIsom[depth-1]:
            for cphase in cphases : 
                G_cand = G_old.copy() 
                G_cand.add_edge(*cphase)
                if not __more_three_multiedges(G_cand): 
                    if not __is_isomorphic_to(G_cand, GNIsom[depth]):
                        GNIsom[depth].append(G_cand)

    # storing stuff 
    LNIG = GNIsom[depth]
    if opt['draw_graphs']+opt['save_edges']: run(['mkdir','-p',opt['dirpath']])

    if opt['save_edges'] : 
        with open(opt['dirpath']+'/edges.json','w') as ouf : 
            json.dump([list(G.edges()) for G in LNIG],ouf)

    if opt['draw_graphs']:
        for i,G in enumerate(LNIG) : 
            GV=pgv.AGraph(directed=False, strict=False)
            GV.add_nodes_from(range(nqubit))
            GV.add_edges_from(list(G.edges()))
            GV.layout()
            GV.draw('%s/graph%i.png'%(opt['dirpath'],i))

    return LNIG



def unique2net(nqubit, net_depth, **kwargs):
    """
    Get a list of 2-bit gates networks. The unique gates are iterated by the following steps:
        1) iterate the non-isomorphic graph (this step doesn't implement parallelism)
        2) place edge ordering from 1). At this step, the relabelling qubit has already implemented.
        3) apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping 

    param
        :nqubit: int, the number of qubits
        :net_depth: int, the depth of the gate-networks
    kwargs
        :path_json: str, path that store a list of edges of non-isomorphic graphs in json
                    format. If this file presents, the iteration of generating new graphs will be started there.
        :draw_graphs: boolean=True, to draw the produced graphs
        :save_edges: boolean=True, save graph as a collection of edges
        :dirpath: str=out-[n]qubit-depth[d], directory path to store outputs

    return:
        a list of 2-bit networks gates, with the LSB convention. For example:
    Example: network (3,5) is 
            q0 -o-o--
                | |
            q1 -o-|--
                  |  
            q2 ---o--
    where (q0,q1)=3, (q0,q2)=5
    """
    #optional arguments
    opt = {'path_json': False, 'draw_graphs': True, 'save_edges': True,
            'dirpath':'out-'+str(nqubit)+'qubit-depth'+str(net_depth)}
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]

    
    start=time()    

    # get a list of non-isomorphic graphs with net_depth edges 
    LNIG = list_non_iso_graphs(nqubit, net_depth, **kwargs)

    # apply ordering 



         

    #print("%i unique network search is obtained within %f seconds"%(len(unique_net),time()-start))

