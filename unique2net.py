#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates criteria of
DiVincenzo and Smolin (cond-mat/9409111). 

The unique gates are iterated by the following steps:
    1) iterate the non-isomorphic graph (non-parallell)
    2) place edge ordering from 1). At this step, the relabelling qubit has already implemented.
    3) apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping 

no parallellization


MAIN USAGE: 

    from unique2net import unique2net

    unique2net(nqubit, network_length)

"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "3.1.2"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"



#standard libraries
from itertools import product, combinations, permutations 
from time import time
from multiprocessing import Process, Value
from numpy import array_split
from math import ceil
from subprocess import run 
from glob import glob
from multiprocessing import Pool, cpu_count


import pygraphviz as pgv
import networkx as nx
import json
import re


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
    """
    Equivalent network by reversing the gates order, basically U^{\dagger}

    :gate_net: tuple(int), a configuration of gate network 

    return tuple
    """
    return gate_net[::-1]   #reverse net


def equiv_bit_permutations(gate_net, nqubit):
    """
    Equivalent networks by bit-relabelling, basically performing identical swaps in
    the beginning and the end, equivalent with simple permutation,
    excluding identity

    :gate_net: tuple(int), a configuration of gate network 
    :nqubit: int, the number of qubits

    return tuple
    """
    result = set([tuple(shuffle_bits(g, p) for g in gate_net) for p in permutations(range(nqubit))])
    result.remove(gate_net) #remove the identity element
    return result

        
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


def equiv_set_net(net, nqubit,**kwargs):
    """
    Get an equivalent set of a network given a network by DS criteria. 

    :net:tuple(int), the gate network
    :nqubit:int, number of qubits

    kwargs : 
        :ds_bit_permutation:boolean=True, include bit permutation criteria 
        :ds_conjugation_by_swap:boolean=True, include conjugation by swapping criteria
        :ds_time_reversal:boolean=True, include time reversal criteria
        :identity:boolean=True, identity is included if True
    """
    #optional arguments
    opt = {'ds_bit_permutation':True, 'ds_conjugation_by_swap':True, 'ds_time_reversal':True, 'identity':True}
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]

    S_equiv = [net] if opt['identity'] else []
    if opt['ds_bit_permutation'] : S_equiv += [*equiv_bit_permutations(net, nqubit)]
    if opt['ds_conjugation_by_swap']: S_equiv += [*equiv_conjugation_by_swapping(net)]
    if opt['ds_time_reversal']: S_equiv += [*equiv_time_reversal(net)]
    
    return {*S_equiv}


## graphs-related helpers
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



def __draw_a_graph(G, outpath='out.png', nodes=False):
    """ Draw a graph

    :G:list(edges) or networkx.MultiGraph or pygraphviy.AGraph
    :outpath:str, the output path
    :nodes:list(int), list of nodes. If it is present, then those nosed will be
           drawn regardless the absence of edges
    """
    GV = pgv.AGraph(directed=False, strict=False)
    if isinstance(G, list):
        GV.add_edges_from(G)
    else:
        GV.add_edges_from(list(G.edges()))

    if nodes : 
        GV.add_nodes_from(nodes)
    GV.layout()
    GV.draw(outpath)


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



def __edges_to_nets_uperm(G, nqubit):
    """
    Return a list of networks from graph up to bit permutation

    :G: networkx.graph or just the edges
    :edges:boolean=False, identifies if the list comprises only edges
    :net:tuple(int), the network configuration
    :nqubit:int, the number of qubits
    """
    net = __edges_to_net(G) 
    GL = list(permutations(net))
    n = len(GL)
    for i in range(n): 
        for j in range(i+1,n):
            equiv_i = equiv_bit_permutations(GL[i],nqubit)
            equiv_j = equiv_bit_permutations(GL[j],nqubit)
            cap = equiv_i.intersection(equiv_j)
            if len(cap):
                GL[i] = False
                break
    return [G for G in GL if G]




def __graphs_to_nets(graph_list, nqubit, edges = False, time_reversal=False):
    """
    Return a set of network from each graph by applying all possible ordering 

    :graph_list:list(nx.MultiGraph) the list of graphs
    :nqubit:int, the number of qubit
    :edges:boolean=False, identifies if the list comprises only edges
    :time_reversal:boolean=False, to include time reversal
    """
    GL = []
    if edges == True : 
            net = __edges_to_net(E)
            GL += __edges_to_nets_uperm(net,nqubit)
    else : 
        for G in graph_list: 
            net = __edges_to_net(G.edges())
            GL +=__edges_to_nets_uperm(net, nqubit)
    GLS = {*GL}
    if time_reversal : 
        GLS.remove(equiv_time_reversal(net))

    return GLS


def __iseqiv_occurrence(net):
    """ Check if the there is more than 3 consecutive identical gates

    :net:tuple(int), the 2-bit gates network 
    """
    count, gate = 0, net[0]  
    for g in net : 
        if g == gate : 
            count += 1
        else : 
            count, gate = 1, g

        if count == 4 : #check >3 conscecutive occurence
            return True

    return False


def draw_equiv_graphs(net, nqubit, criteria ,images_per_row=4, dirpath='out'):
    """
    Get a pciture comprises a set of equivalent networks by a criteria of DS 

    :net:tuple(int), the network
    :nqubit:int, the number of qubits
    :criteria: 'permutation_bit' | 'conjugation_by_swapping' | 'time_reversal'
    :images_per_row:int, the number of graph displayed per row
    :outpath:str, the path for the output
    """
    if criteria[0] == 'p':
        equivs, idf = equiv_bit_permutations(net, nqubit), 'equiv_perm'
    elif criteria[0] == 'c':
        equivs, idf = equiv_conjugation_by_swapping(net), 'equiv_conj'
    else:
        equivs, idf = {equiv_time_reversal(net)}, 'equiv_trev'
    
    netidf = re.sub('\ |\(|\)','',str(net))#network identifier
    netidf = re.sub(',','-',netidf)

    run(['mkdir','-p','%s/%s/%s'%(dirpath, idf, netidf)])

    __draw_a_graph(__net_to_graph(net), '%s/%s/%s/00.png'%(dirpath,idf,netidf), range(nqubit))
    for i, n in enumerate(equivs):
        __draw_a_graph(__net_to_graph(n), '%s/%s/%s/%i.png'%(dirpath,idf,netidf,i), range(nqubit))

    # combine graphs per row, then combine all rows
    fnames = array_split(glob('%s/%s/%s/*.png'%(dirpath,idf,netidf)),max(ceil(len(equivs)/images_per_row),1))
    for i,row in enumerate(fnames):
        run(['convert', *list(row), '+append', '-border','20','%s/%s/%s/row%i.png'%(dirpath,idf,netidf,i)])

    # combine all graph and remove the image per graph 
    run(['convert', *glob('%s/%s/%s/row*.png'%(dirpath,idf,netidf)),'-append','-border','20','%s/%s/%s.png'%(dirpath,idf,netidf)])
    run(['rm', '-r',  '%s/%s/%s'%(dirpath,idf,netidf)])




def __worker_draw_equivs(net, nqubit, dirpath, images_per_row):
    """
    worker that draw equivalent nodes based on criterias DS: bit permutation and conjugation by swapping 

    :net:tuple(int), the network
    :nqubit:int, the number of qubit
    :outpath:str, the path directory
    :images_per_row:int, the number of graph displayed per row
    """
    run(['mkdir','-p',dirpath])
    draw_equiv_graphs(net,nqubit, 'p' ,images_per_row=images_per_row, dirpath=dirpath)
    draw_equiv_graphs(net,nqubit, 'c' ,images_per_row=images_per_row, dirpath=dirpath)
    draw_equiv_graphs(net,nqubit, 't' ,images_per_row=images_per_row, dirpath=dirpath)



def __worker_draw_equiv_mapper(args):
    """ multi-arguments helper for __worker_draw_equivs
    """
    return __worker_draw_equivs(*args)



### test methods

def test_draw_by_equivalents(nqubit, net_depth, ncpu=4, images_per_row=4, dirpath='out'): 
    """
    Test function: draw graph by grouping its equivalents. Each folder contains
    equivalent graphs grouped by DS equivalents

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks
    :ncpu:int, number of cpu for parallelization
    :images_per_row:int, the number of graphs drawn per row
    :dirpath:str="out", output path 
    """
    run(['mkdir','-p',dirpath]) 
    cphases = [i for i in range(2**nqubit) if bin(i).count('1')==2]
    Lstart = [net for net in product(cphases, repeat=net_depth) if not __iseqiv_occurrence(net)]
    Largs = [(net, nqubit, dirpath ,images_per_row) for net in Lstart]
    P = Pool(ncpu)
    P.map(__worker_draw_equiv_mapper, Largs)
                


def is_equiv_iso_conj(net, non_iso_graphs):
    """
    Tell if a network is equivalent by criteria "conjugation by swapping" using graph isomorphism.

    :net:tuple(int), the network
    :non_iso_graphs:iterative(networkx.MultiGraph), a list of non isomorphic graphs
    """
    G_test, answers = __net_to_graph(net, 'multigraph'), 0

    for G in non_iso_graphs:
        answers += nx.is_isomorphic(G,G_test) 
        if answers > 2 : return True

    if answers > 2: 
        return True
    else : 
        return False



### main methods

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
            __draw_a_graph(G, '%s/graph%i.png'%(opt['dirpath'],i))

    return LNIG


def __helper_edges_to_nets_uperm(args):
    return __edges_to_nets_uperm(*args)

def __release_list(L):
    del L[:]
    del L

def __net_equiv_conj(net, i):
    """networks that have equivalence by conjugation by swapping 
    :net: tuple, the network gates
    :i:int, the identifier

    return dict {index, set_of_equivalent_nets}
    """
    nets_conj = dict()  
    netsc = equiv_conjugation_by_swapping(net)
    if len(netsc) > 0 : 
        nets_conj[i] = netsc
        return nets_conj
    else : return False

def __helper_net_equiv_conj(args):
    return __net_equiv_conj(*args)


def __toel_conj(conj_dnets, LNIG):
    """ to eliminate ones that are equivalent by conjugation by swap

    :conj_equiv_nets: equivalent networks
    :LNIG:list of non isomorphic graphs
    """
    nets = list(conj_dnets.values())[0]
    for net in nets: 
        if is_equiv_iso_conj(net, LNIG):
            return list(conj_dnets.keys())[0]
    return False


def __helper_toel_conj(args):
    return __toel_conj(*args)


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
        :ncpu:int=cpu_count, the number of cpu in parallelization 
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
            'ds_bit_permutation': False, 'ds_conjugation_by_swap':True,
            'ds_time_reversal': False, 'ncpu':cpu_count() }
    for key in opt: 
        if key in kwargs : opt[key]=kwargs[key]

    start=time()    

    # get a list of non-isomorphic graphs with net_depth edges 
    list_niso_kwargs = {k: opt[k] for k in ['path_json','draw_graphs','dirpath']} 
    LNIG = list_non_iso_graphs(nqubit, net_depth, **list_niso_kwargs)

    # apply ordering, we obtain unique networks up to permutation 
    P = Pool(opt['ncpu'])
    LNets = P.map(__helper_edges_to_nets_uperm,[(G.edges(), nqubit) for G in LNIG])
    Nets = [] #flattened
    for N in LNets : Nets.extend(N)
    __release_list(LNets)

    #get ones that have equivalences in conjugation by swap
    P = Pool(opt['ncpu'])
    nets_conj = P.map(__helper_net_equiv_conj,[(net,i) for i,net in enumerate(Nets)] ) 
    nets_conj = filter(lambda x:x, nets_conj)

    #then, eliminate
    P = Pool(opt['ncpu'])
    toel = P.map(__helper_toel_conj, [(dnet, LNIG) for dnet in nets_conj])

    for i in toel : 
        if i : Nets[i] = False
                
    Nets=[net for net in Nets if net]                
    print("%i unique network search is obtained within %f seconds"%(len(Nets),time()-start))

    return Nets



