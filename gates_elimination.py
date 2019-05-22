#!/usr/bin/env python3
from itertools import product, combinations #standard lib

""" gates_elimination.py: list all unique gates by eliminations
    
    arXiv ref: cond-mat/9409111.
"""

__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"




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

    return int
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



## listing gates and networks ##

def all_2g(nqubit):
    """
    Return a set of all possible 2-qubit gates of nqubits qubits. Those
    will be a set of integers composing the network gates

    :nqubit: int, the number of qubits 

    return set
    """
    return set([a for a in range(2**nqubit) if bin(a).count('1')==2])



def all_2g_networks(netsize, gates2):
    """
    Return a list composing tuples of all possible 2-bit gates network
    with size netsize

    :nqubit: int, the number of qubits 
    :gates2: set, the set 2-qubit gates compsing the network

    return list
    """
    return list(product(gates2, repeat=netsize))


## eliminations of 2-bit gate networks ##

def eliminate3(gate_list):
    """
    Eliminate the networks with > 3 consecutive occurrence of CNOTS. 
    Example: (1,3,3,3,3,1) is eliminated

    :gate_list: list(tuple), the list of gate networks 

    return list
    """
    for i, net in enumerate(gate_list) : 
        count, gate = 0, net[0]  
        for g in net : 
            if g == gate : 
                count += 1
            else : 
                gate = g
                count = 1

            if count > 3 : #check >3 conscecutive occurence
                gate_list[i] = False
                break

    #filter out the False nets
    return [g for g in gate_list if g]



def eliminate_time_reversal(gate_list):
    """
    Eliminate the ones in the reversed order

    :gate_list: list(tuple), the list of gate networks

    return list
    """
    for i, net in enumerate(gate_list):
        if net[::-1] in gate_list : 
            gate_list[i] = False

    #filter out the False nets
    return [g for g in gate_list if g]



def eliminate_relabelling(nqubit, gate_list):
    """
    Elimiate the ones which are equivalent by bit-relabelling
    It is equivalent by performing one swap

    :nqubit: int, the number of qubits 
    :gate_list: list(tuple), the list of gate networks

    return list(tuple)
    """
    #all possible swaps 
    lswaps = combinations(range(nqubit), 2) 

    for swap in lswaps:
        for i, net in enumerate(gate_list):
            net2 = []   #the swapped net
            if net : 
                for g in net : #construct the new net 
                    g2 = swapbits(*swap, g)   
                    net2 += [g2]
                if tuple(net2) in gate_list : 
                    gate_list[i] = False

    #filter out the False nets
    return [g for g in gate_list if g]
        
        
def eliminate_conjugation_by_swapping(gate_list):
    """
    Eliminate conjugation by swapping. Choose the a gates sandwiched by
    identical gates, then swap it.

    :gate_list: list(tuple), the list of gate networks

    return list(tuple)
    """
    for i, net in enumerate(gate_list):
        for j, g in enumerate(net[1:-1]) : #inside big sandwich
            mnet = list(net)  #get a mutable object
            g1, g2 = net[j], net[j+2]  #small sandwich g1|g|g2
            if g1 == g2 :
                poss = get_pos_ones(g1)
                g_swapped = swapbits(*poss, g) 
                mnet[j+1] = g_swapped

                if tuple(mnet) in gate_list : 
                    gate_list[i] = False

    #filter out the False nets
    return [g for g in gate_list if g]



## main function ##

def get_unique_2gates_networks(nqubit, network_length):
    """
    Get a list of unique 2-bit gates network by four eliminations
    The format is binary with LSB convention. 

    example: 
            q0 ----
            q1 ---- 
            q2 ----
    gate (q0,q1)=3, (q0,q2)=5, (q1,q2)=6 

    :nqubit: int, the number of qubits
    """
    g2 = all_2g(nqubit) #all kind of gates

    L = all_2g_networks(network_length, g2)
    print(len(L))
    L = eliminate3(L)
    L = eliminate_relabelling(nqubit, L)
    L = eliminate_time_reversal(L)
    L = eliminate_conjugation_by_swapping(L)
    print(len(L))

    print("second iteration")

    L = eliminate3(L)
    L = eliminate_relabelling(nqubit, L)
    L = eliminate_time_reversal(L)
    L = eliminate_conjugation_by_swapping(L)
    print(len(L))

    
    return L



