#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates by eliminations from
DiVincenzo and Smolin (cond-mat/9409111). 
This version 2 applies different approach. Instead of eliminating the
equivalent ones, the unique ones are added one by one.

Parallelization is implemented.


    
    MAIN USAGE: 

        from unique2net import unique2net

        unique2net(nqubit, network_length, NCPU)

"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "2.1.0"
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




## listing gates and networks ##

def all_2g(nqubit):
    """
    Return a set of all possible 2-qubit gates of nqubits qubits. Those
    will be a set of integers composing the network gates

    :nqubit: int, the number of qubits 

    return set
    """
    return set([a for a in range(2**nqubit) if bin(a).count('1')==2])



def net_eliminate3(gate_net):
    """
    Eliminate the networks with > 3 consecutive occurrence of CNOTS. 
    Example: (1,3,3,3,3,1) is eliminated

    :gate_net: tuple(int), a configuration of gate network 

    return boolean
    """
    count, gate_counted = 0, gate_net[0]  
    for g in gate_net : 
        if g == gate_counted : 
            count += 1
        else : 
            gate_counted = g
            count = 1

        if count == 4 : #check >3 conscecutive occurence
            return False
    return True


def all_2g_networks(net_depth, gates2):
    """
    Return a list composing tuples of all possible 2-bit gates(CPHASEs) network
    with size netsize. The returned networks has no more than
    consecutive occurance more than 3 CPHASE.

    :net_depth: int, the network depth, basically the number of CPHASE 
    :gates2: set, the set 2-qubit gates compsing the network

    return filter
    """
    return filter(net_eliminate3, product(gates2, repeat=net_depth))


## obtaining equivalent networks network ##


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


def equiv_DS(nqubit, gate_net):
    """
    Get a set of all equivalent networks based on
    DiVincenzo and Smolin equivalent networks

    :nqubit: int, the number of qubits used
    :gate_list: list(tuple|bool), the list of gate networks
    """
    equivs = set()
    #equivs = equivs.union([equiv_time_reversal(gate_net)])
    equivs = equivs.union(equiv_conjugation_by_swapping(gate_net))
    equivs = equivs.union(equiv_bit_permutations(gate_net, nqubit))
    
    return equivs


## formating 

def change_format(g2,nqubit,G):
    """
    This function changes the format from binary type to tuple type
    to the index in combo
    
    Eg. If the gate number in binary is 6=110, the tuple comes out as (2,3)
    (i.e. get_pos_one(6) + 1) which further is changed to the index of that tuple
    in c i.e 2. 
    """
    a=list(range(1,nqubit+1))
    c=list(combinations(a,2))
    g2_dict=dict()
    
    for i in list(g2):
        pos_tuple = tuple(x+1 for x in get_pos_ones(i))
        g2_dict[i]=c.index(pos_tuple)
        
    initial_list =list(G)
    final_list=[]
    for i in initial_list:
        new_tuple = tuple(g2_dict[x] for x in tuple(i))
        final_list.append(new_tuple)
    
    return(final_list)



## main function ##
def __worker_net_checker(nqubit, unet_list, net_test, found):
    """
    Worker to check if the given gate network is in the list

    :nqubit: int, the number of qubits
    :unet_list: list(tuple), the list of unique gates that will be compared to
                the net_test  
    :net_test: tuple, the tested gate network
    :found: multiprocessing.Value ctype int, the value  to track if the equivalent
            network is found
    """ 
    net_test_equiv = equiv_DS(nqubit, net_test)
    for unet in unet_list : 
        if found.value:
            break
        if equiv_DS(nqubit, unet).intersection(net_test_equiv) : #if equivalent
            with found.get_lock():
                found.value += 1
            break



def elimination_methode2():
    """
    Get a list of unique 2-bit gates network by four eliminations
    The format is binary with LSB convention. 
    Example: 
            q0 ----
            q1 ---- 
            q2 ----
    gate (q0,q1)=3, (q0,q2)=5, (q1,q2)=6 

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks
    :NCPU: int, the number of cpu in multiprocessing

    return list
    """
    unique_net = []
    for net in all_2g_networks(net_depth, all_2g(nqubit)): 

        unet_idx = array_split(range(len(unique_net)), NCPU) #to split the unique nets wrt indices 
        proc_num = min(NCPU, len(unique_net))

        found = Value('i', 0)
        if proc_num : 
            found.value = 0
            pool = [Process(target=__worker_net_checker, 
                args=(nqubit,
                    unique_net[unet_idx[i][0]:unet_idx[i][-1]+1], net,
                    found)) for i in range(proc_num)]

            ## parallel section
            for procc in pool : 
                procc.start() #spawn each thread
            for p in pool: #join each thread 
                p.join()
            ## end of parallel section
        
        if not found.value: 
            unique_net.append(net)


    return unique_net


def is_isomorphic_to(G_test, G_list):
    """
    Check if G_test is isomorphic to one of the graph in G_list

    :G_test: nx.graph, the graph to be tested
    :G_list: list(nx.graph), the list of graph to be compared to
    """
    for G in G_list : 
        if nx.is_isomorphic(G,G_test):
            return True
    return False



def ListNonIsoGraphs(nqubit, net_depth, *path_dict):
    """
    List the non-isomorphic graphs

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks
    :path_dict: str, path that store a list of non-isographic graphs 
    """
    cphases = list(combinations(range(nqubit),2)) 

    #Non-Isomorphic graphs
    GNIsom = dict()
    GNIsom[0]=[nx.MultiGraph()]
    GNIsom[0][0].add_nodes_from(range(nqubit))
    GNIsom[1]=[nx.MultiGraph()]
    GNIsom[1][0].add_edge(*cphases[0])

    for depth in range(2, net_depth+1):
        GNIsom[depth]=list()
        for G_old in GNIsom[depth-1]:
            for cphase in cphases : 
                G_cand = G_old.copy() 
                G_cand.add_edge(*cphase)
                if not is_isomorphic_to(G_cand, GNIsom[depth]):
                    GNIsom[depth].append(G_cand)
    return GNIsom[net_depth]




def graph_methode(nqubit, net_depth):
    """
    Get a list of 2-bit gates network using graph isomorphism

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks
    """
    outpath='out-%iqubit-depth%i'%(nqubit,net_depth)
    run(['mkdir','-p',outpath])
    LNIG = ListNonIsoGraphs(nqubit, net_depth)

    for i,G in enumerate(LNIG) : 
        GV=pgv.AGraph(directed=False, strict=False)
        GV.add_edges_from(list(G.edges()))
        GV.layout()
        GV.draw('%s/graph%i.png'%(outpath,i))
         




def unique2net(nqubit, net_depth, methode=graph_methode):
    """
    Get a list of unique 2-bit gates network by four eliminations
    The format is binary with LSB convention. 
    Example: 
            q0 ----
            q1 ---- 
            q2 ----
    gate (q0,q1)=3, (q0,q2)=5, (q1,q2)=6 

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks
    :methode: function, the method to do find the unique networks
    :method_kwargs: dict, arguments for the methode

    return list
    """
    start=time()    
    unique_net = methode(nqubit, net_depth)
    #print("%i unique network search is obtained within %f seconds"%(len(unique_net),time()-start))




