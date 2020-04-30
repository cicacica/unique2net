#!/usr/bin/env python3

__doc__=""" Contains classes:

GraphQGate: creates object that contains graph including its gate ordering. This
object is to be associated with graph and quantum gates.

Net2Gates: creates an object for the 2-bit gates network

bitop: class contains staticmethods for bit-operations related
"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"



#standard libraries
import pygraphviz as pgv
import networkx as nx
import os
import re

from subprocess import run
from multiprocessing import Pool, cpu_count
from glob import glob
from numpy import array_split
from math import ceil, log10
from itertools import combinations, groupby
from more_itertools import consecutive_groups



# Due to restriction of multiprocessing
def _GraphQNet__helper_draw_a_graph(args):
    return _GraphQNet__draw_a_graph(*args)

def _GraphQNet__draw_a_graph(nqubit, netgates, outdir, fname):
    ng = GraphQNet(nqubit,netgates)
    ng.set_out_dir(outdir)
    ng.draw_graph(fname)


class GraphQNet:
    """
    Network of 2-qubit gates object
    """
    def __init__(self, nqubit, netgates):
        """ Instantiation of the GraphQNet object. That is the object contains
        a 2-bit gates network and can has some graph representation.

        :nqubit: int, the number of qubits.
        :netgates: tuple(int), the gate network
        """
        # input checks
        if type(nqubit) != int :
            raise TypeError('nqubit -- number of qubits -- requires int')
        if type(netgates) == tuple :
            if type(netgates[0]) == int :
                pass
        else :
            raise TypeError('netgates is a tuple of integers')

        self.nqubit = nqubit
        self.netgates = netgates
        self.outdir = False
        self.graph = False
        self.depth = len(netgates)

        #set some stuff
        self.set_out_dir()
        self.set_graph()

    def __copy__(self):
        return GraphQNet(self.nqubit, self.netgates)

    @staticmethod
    def edges_to_net(edges):
        """ Return conversion of a graph to a bit network convention with an arbitrary
            2-bit gates ordering

            :edges: list(tuple), the list of edges
        """
        return tuple([bitop.pos_ones_toint(*e) for e in edges])

    def add_ordered_edge(self, edge):
        """ Add an edge with ordering as number of edge --- ordering starts from 0

        :edge:tupe(int,int) the edge
        """
        self.graph.add_edge((*edge, self.depth), weight='ordering')
        self.depth += 1
        self.netgates = (*self.netgates, bitop.pos_ones_toint(*edge))

    @staticmethod
    def net_to_edges(netgates):
        """ Return the conversion of a 2-bit gate network to edges with nodes 0 ... n

            :netgates: tuple(int), two-bit network netgates
        """
        return [bitop.pos_of_ones(gate) for gate in netgates]


    def set_netgates(self, new_netgates):
        """ Renew the netgates attribute
        """
        self.netgates = new_netgates
        self.depth = len(new_netgates)
        self.set_graph()

    def set_graph(self):
        """ Set self.graph
        """
        wedges = [(*e, i) for i,e in enumerate(self.net_to_edges(self.netgates))]
        self.graph = nx.MultiGraph()
        self.graph.add_weighted_edges_from(wedges, weight='ordering')
        if len(self.graph.nodes()) > self.nqubit :
            raise ValueError('Hi there, you need at least %i qubits'%len(self.graph.nodes()))
        self.graph.add_nodes_from(range(self.nqubit))


    def set_out_dir(self, *outdir):
        """
        Set the output directory

        :outdir: str, the output directory
        """
        self.outdir = outdir[0] if outdir else os.getcwd()
        run(['mkdir', '-p', self.outdir])


    def more_three_con_edges(self):
        """ Tells if the network has more than three consecutive edges
        """
        for key, l in groupby(self.netgates):
            if len(list(l)) > 3 :
                return True
        return False

    def conjugation_by_swap(self):
        """ Return a set of GraphQNet objects, the equivalent networks by swap conjugation, if there is any.
        """
        #group the element by occurences, only the ones that
        #occur more than once, has potential to be a sanwdich
        sw = [k for k,l in groupby(sorted(self.netgates)) if len(list(l))>1]

        #a list of contains network, and index where it starts
        sidx, i = [], 0
        for k, l in groupby(self.netgates):
            if k in sw:
                sidx += [(k, i)]
            i += len(list(l))

        #start index by keys, sorted already
        sikey = dict([(k,[]) for k in sw ])
        for k, idx in sidx :
            sikey[k].append(idx)

        #get all combinations of sandwiches
        #unique up to the conjugation by swap
        unet = []
        for k, sidx in sikey.items():
            idxs =[(a,b) for a,b in combinations(sidx,2) if (b-a)>1]
            s1,s2 = bitop.pos_of_ones(k)
            for i1, i2 in idxs :
                #swap indices between i1 xxx i2, by sandwiching with swap k
                new_net = list(self.netgates)
                for i in range(i1+1, i2):
                    new_net[i] = bitop.swap(self.netgates[i],s1,s2)
                unet.append(tuple(new_net))

        unet = set(unet)
        if self.netgates in unet :
            unet.remove(self.netgates)

        gqn_list = [GraphQNet(self.nqubit, net) for net in unet]

        return gqn_list


    def time_reversal(self):
        """
        Return network, which is it's time-reversal
        equivalence: equivalent network by reversing the gates order, basically
        U^{\dagger}
        """
        return GraphQNet(self.nqubit, self.netgates[::-1])

    def draw_graph(self, outfile='graphnet.png'):
        """
        :outfile:output file
        """
        gv = pgv.AGraph(directed=False, strict=False)

        for node in range(self.nqubit):
            gv.add_node(node, shape='circle')
        for i,gate in enumerate(self.netgates) :
            gv.add_edge(*bitop.pos_of_ones(gate), label=str(i))
        gv.layout()
        gv.draw('%s/%s'%(self.outdir, outfile))


    @staticmethod
    def __compare_edges(G1, G2):
        return G1[0]['ordering'] == G2[0]['ordering']


    def is_isomorphic_to(self, GQN):
        """
        Check if G_test is isomorphic to another GraphQNet instance.
        It includes bit-permutation and conjugation by swap in DS criteria
        """
        return  nx.is_isomorphic(self.graph, GQN.graph, edge_match=self.__compare_edges)

    def is_isomorphic_uptolist(self, list_gqn):
        """
        Check if the network is isomorphic, compared to every
        element in list_gqn.
        :list_gqn:list(GraphQNet) list of objects GraphQNet
        """
        for gqn in list_gqn:
            if self.is_isomorphic_to(gqn):
                return True
        return False


    @staticmethod
    def draw_netgraphs_list(netgates_list, nqubit, images_per_row=False, outfile='picture.png', nrow=False, ncpu=False):
        """
        Get a picture contains graphs, where each graph is in netgates_list.

        :netgates_list:list(tuple(int))  the gate networks
        :nqubit:int, the number of qubits
        :images_per_row:False, the number of graph displayed per row
        :outpath:str, the path for the output
        :nrow:int, the number of row of image tiles
        :ncpu:int, the cpu number
        """
        def get_netid(net):
            #get network identifier as filenames
            netid = re.sub('\ |\(|\)','',str(net))#network identifier
            return re.sub(',','-',netid)

        outdir = os.getcwd()+'/'+str(os.getpid())
        run(['mkdir', '-p', outdir])
        ncpu = ncpu if ncpu else cpu_count()

        fnames = [get_netid(net)+'.png' for net in netgates_list]
        args = [(nqubit, net, outdir, fname) for net,fname in zip(netgates_list,fnames)]

        P = Pool(ncpu)
        P.map(__helper_draw_a_graph, args)
        P.close()
        P.join()

        fpaths = glob(outdir+'/*.png')
        num_row = nrow if nrow else max(int(ceil(len(args)/8)),1) 
        run(['montage', *list(fpaths), '-mode', 'Concatenate', '-tile' ,  'x%i'%num_row,  outfile ])

   #    # combine all graphs per row, then combine all rows
   #    images_per_row = images_per_row if images_per_row else max(int(log10(len(args))),1)*5
   #    fpathss = array_split(glob(outdir+'/*.png') ,max(ceil(len(args)/images_per_row),1))
   #    for i,row in enumerate(fpathss):
   #        run(['convert', *list(row), '+append', '-border','2','-alpha','set','%s/row%i.png'%(outdir,i)])

   #    # combine all graph and remove the image per graph
   #    run(['convert', *glob('%s/row*.png'%outdir),'-append','-border','2', '-alpha','set',outfile])
        run(['rm', '-r', outdir])


class bitop:
    """
    Bit-operation related methods
    """

    @staticmethod
    def bit_at(num, k):
        """
        get the k-th bit of num, where num has length nqubit

        :num: int, the identified number
        :k: int, the position

        return 0 or 1
        """
        return (num & (1 << k)) >> k

    @classmethod
    def swap(cls, num, p1, p2):
        """
        Swap p1-th bit with p2-th bit within number num

        :num: int, the number
        :p1:,:p2: int, the swapped position

        return int"
        Shuffle the binary of num with some permutation

        :num: int, the number
        :permutaion: tup
        """
        b1 = cls.bit_at(num, p1)
        b2 = cls.bit_at(num, p2)

        # XOR the two
        xor = b1^b2

        # get a mask xor with position 00xor000xor00
        #                                p1    p2
        xor = (xor << p1) | (xor << p2)

        # or swap with the original number
        # if xor=0, number is the same, otherwise 0->1 1->0
        res = num ^ xor

        return res


    @staticmethod
    def pos_of_ones(num):
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


    @classmethod
    def shuffle(cls, num, permutation):
        """
        Shuffle the binary of num with some position permutation

        :num: int, the number
        :permutation: tuple(int), the final permutation
        """
        num2 = 0
        for i, new_pos in enumerate(permutation):
            # updating from right to left 0xxx | j000
            num2 =  num2 | (cls.bit_at(num, new_pos) << i)

        return num2


    @staticmethod
    def pos_ones_toint(*args):
        """
        Return the integer given posisition of ones of a bitstring

        :args:int, the posisions of ones
        """
        res = 0
        for a in args :
            res += 2**a
        return res

