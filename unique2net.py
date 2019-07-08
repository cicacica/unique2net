#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates by eliminations from
DiVincenzo and Smolin (cond-mat/9409111). 
This version 2 applies different approach. Instead of eliminating the
equivalent ones, the unique ones are added one by one.

Parallelization is not yet implemented.


    
    MAIN USAGE: 

        from unique2net import unique2net

        unique2net(nqubit, network_length)

"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "2.0.0"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"



#standard libraries
from itertools import product, combinations, permutations 
from time import time




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


def equiv_bit_permutations(gate_net):
    """
    Equivalent networks by bit-relabelling, basically performing identical swaps in
    the beginning and the end, equivalent with simple permutation,
    including identity

    :gate_net: tuple(int), a configuration of gate network 

    return tuple
    """
    return [tuple(shuffle_bits(g, p) for g in gate_net) 
            for p in permutations(range(len(gate_net)))]

        
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


def equiv_DS(gate_net):
    """
    Get all equivalent networks based on
    DiVincenzo and Smolin equivalent networks

    :nqubit: int, the number of qubits used
    :gate_list: list(tuple|bool), the list of gate networks
    """
    equivs = set()
    equivs = equivs.union([equiv_time_reversal(gate_net)])
    equivs = equivs.union(equiv_conjugation_by_swapping(gate_net))
    equivs = equivs.union(equiv_bit_permutations(gate_net))
    
    return equivs



## main function ##

def __compare_net_to_equiv(net, equivalents):
    """
    Helper of unique2net function. It compares a network among given a
    list of equivalent networks. The equivalent networks of net is
    listed using DiVincenzo - Smolin equivalent
    criteria 

    :net: tuple(int), a configuration of gate network 
    :equivalents: list(tuple), a list of gate network 

    return boolean
    """
    for n in equiv_DS(net) : 
        if n in equivalents : 
            return True
    return False
    


def unique2net(nqubit, net_depth):
    """
    Get a list of unique 2-bit gates network by four eliminations
    The format is binary with LSB convention. 

    example: 
            q0 ----
            q1 ---- 
            q2 ----
    gate (q0,q1)=3, (q0,q2)=5, (q1,q2)=6 

    :nqubit: int, the number of qubits
    :net_depth: int, the depth of the gate-networks

    return generator
    """
    unique_net = []

    start=time()    
    for net in all_2g_networks(net_depth, all_2g(nqubit)): 
        equiv1 = equiv_DS(net)
        exists = False
        for unet in unique_net : 
            if __compare_net_to_equiv(unet, equiv1) : 
                exists = True
                break
        if not exists : 
            unique_net.append(net)

    print("%i unique network search is obtained within %f seconds"%(len(unique_net),time()-start))

    return unique_net





