#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates criteria of
DiVincenzo and Smolin (cond-mat/9409111). 

The unique gates are iterated by the following steps:
    1) iterate the non-isomorphic graph 
    2) place edge ordering from 1). At this step, the relabelling qubit has already implemented.
    3) apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping 

no parallellization


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
    """Equivalent network by reversing the gates order, basically U^{\dagger}

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
    return set([tuple(shuffle_bits(g, p) for g in gate_net) 
            for p in permutations(range(nqubits))])

        
def equiv_conjugation_by_swapping(gate_net):
    """
    Equivalent network that is equivalent by swapping conjugation.
    Choose the a gates sandwiched by identical gates, then swap it.

    :gate_net: tuple(int), a configuration of gate network 

    return set(tuple)
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
                result += [tuple(mnet)]

    return set(result)


def nets_with_ds_equivalents(nets, nqubit,**kwargs):
    """
    convert a list of networks into a list of classes of equivalent networks.

    :nets: list(tuple), the list of networks
    :nqubit:int, the number of nqubits
    kwargs
        :ds_bit_permutation:boolean=False, include bit permutation criteria 
        :ds_time_reversal:boolean=False, include time reversal criteria
        :ds_conjugation_by_swap:boolean=True, include conjugation by swapping criteria
    """
    opt = {'ds_bit_permutation': False, 'ds_conjugation_by_swap':True, 'ds_time_reversal': False }
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]

    nets_equivs=list()
    if opt['ds_bit_permutation'] and not opt['ds_conjugation_by_swap'] and not opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_bit_permutations(net,nqubit)]
    elif not opt['ds_bit_permutation'] and opt['ds_conjugation_by_swap'] and not opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_conjugation_by_swapping(net)]
    elif not opt['ds_bit_permutation'] and not opt['ds_conjugation_by_swap'] and opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_time_reversal(net)]
    elif opt['ds_bit_permutation'] and opt['ds_conjugation_by_swap'] and not opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_bit_permutations(net,nqubit).union(equiv_conjugation_by_swapping(net))]
    elif opt['ds_bit_permutation'] and not opt['ds_conjugation_by_swap'] and opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_bit_permutations(net,nqubit).union(equiv_time_reversal(net))]
    elif opt['ds_bit_permutation'] and opt['ds_conjugation_by_swap'] and not opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_bit_permutations(net,nqubit).union(equiv_conjugation_by_swapping(net))]
    elif not opt['ds_bit_permutation'] and opt['ds_conjugation_by_swap'] and opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_conjugation_by_swapping(net).union(equiv_time_reversal(net))]
    elif opt['ds_bit_permutation'] and not opt['ds_conjugation_by_swap'] and opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_bit_permutations(net,nqubit).union(equiv_time_reversal(net))]
    elif not opt['ds_bit_permutation'] and opt['ds_conjugation_by_swap'] and opt['ds_time_reversal']:
        for net in nets: nets_equivs += [equiv_conjugation_by_swapping(net).union(equiv_time_reversal(net))]
    else :
        for net in nets : nets_equivs += [equiv_bit_permutations(net,nqubit).union(quiv_conjugation_by_swapping(net)).union(
                                        equiv_time_reversal(net))]
    return nets_equivs


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
    List the non-isomorphic graphs or edges, where no more that three edges between two nodes.

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
    opt = {'path_json': False, 'draw_graphs': True, 'save_edges': True,
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

        if start_depth > net_depth : 
            return list_edges #nothing to do here  

        GNIsom[start_depth-1] = [nx.MultiGraph(edges) for edges in list_edges]
        print("iteration starts from %i edges..."%(start_depth-1))

    else : 
        GNIsom[1]=[nx.MultiGraph([cphases[0]])]
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



# listing nets related methods
def __edges_to_net(edges):
    """
    convert edges of a graph to a bit network with an arbitrary 2-bit gates ordering

    :edges: list(tuple), the list of edges
    """
    return tuple([2**a+2**b for a,b in edges])
    

def __net_to_graph(network, graph_type='graphviz'):
    """
    convert a network into a graph
    
    :network:tuple, the 2-qubit network gates with lsb convention
    :graph_type:str, 'graphviz' or 'multigraph'
    """
    if graph_type[0] == 'm':
        return nx.MultiGraph([get_pos_ones(gate) for gate in network]) 
    else : 
        G = pgv.AGraph(directed=False, strict=False)
        G.add_edges_from([get_pos_ones(gate) for gate in network]) 
        return G


def __graphs_to_nets(graph_list, edges = False):
    """
    Return a list of network from each graph by applying all possible ordering 

    :graph_list:list(nx.MultiGraph) the list of graphs
    :edges:boolean=False, identifies if the list comprises only edges
    """
    GL = []

    if edges == True : 
        for E in graph_list: 
            net = __edges_to_net(E)
            GL += list(set(permutations(net)))
    else : 
        for G in graph_list: 
            net = __edges_to_net(G.edges())
            GL += list(set(permutations(net)))

    return GL



def unique2net(nqubit, net_depth, **kwargs):
    """
    Get a list of 2-bit gates networks. The unique gates are iterated by the following steps:
        1) iterate the non-isomorphic graph (this step doesn't implement parallelism)
        2) place edge ordering from 1). At this step, the relabelling qubit has already implemented.
        3) apply more criteria of DiVincenzo and Smolin: time reversal and conjugation by swapping 

    param
        :nqubit: int, the number of qubits
        :net_depth: int, the depth of the gate-networks
    kwargs
        :path_json: str, path that store a list of edges of non-isomorphic graphs in json
                    format. If this file presents, the iteration of generating new graphs will be started there.
        :draw_graphs: boolean=True, to draw the produced graphs
        :save_edges: boolean=True, save graph as a collection of edges
        :dirpath: str=out-[n]qubit-depth[d], directory path to store outputs
        :ds_bit_permutation:boolean=False, include bit permutation criteria 
        :ds_time_reversal:boolean=False, include time reversal criteria
        :ds_conjugation_by_swap:boolean=True, include conjugation by swapping criteria

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
            'dirpath':'out-'+str(nqubit)+'qubit-depth'+str(net_depth),
            'ds_bit_permutation': True, 'ds_conjugation_by_swap':True,
            'ds_time_reversal': False }
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]

    start=time()    

    # get a list of non-isomorphic graphs with net_depth edges 
    list_niso_kwargs = {k: opt[k] for k in ['path_json','draw_graphs','dirpath']} 
    LNIG = list_non_iso_graphs(nqubit, net_depth, **list_niso_kwargs)

    # apply ordering 
    LNets = __graphs_to_nets(LNIG, type(LNIG[0]==list))

    # get nets with it's equivalents
    ds_kwargs = {k:opt[k] for k in ['ds_bit_permutation','ds_time_reversal','ds_conjugation_by_swap']}
    net_w_equiv = nets_with_ds_equivalents(LNets, nqubit,**ds_kwargs)
    
    for i, s_equiv in enumerate(net_w_equiv) : 
        if len(s_equiv): 
            for k,net in enumerate(LNets) :  
                if i!=k and net in s_equiv : 
                    LNets[i]=False
                    break
    unique_net = [net for net in LNets if net]
    print("%i unique network search is obtained within %f seconds"%(len(unique_net),time()-start))

    return unique_net



