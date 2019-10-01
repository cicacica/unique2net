#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates criteria of
DiVincenzo and Smolin (cond-mat/9409111).

The unique gates are iterated by the following steps:
    1) iterate the non-isomorphic graph (non-parallell)
    2) place edge ordering from 1). At this step, the relabelling qubit has already implemented.
    3) apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping


MAIN USAGE:

    from unique2net import unique2net

    unique2net(nqubit, network_length)

"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "4.0.0"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"



#standard libraries
from itertools import combinations
from time import time
from subprocess import run
from multiprocessing import cpu_count
import json

#additional library
from GraphQNet import GraphQNet, bitop




def iterate_graphqnet_noniso(nqubit, graphqnet_list, net_edges):
    """
    Produce non-isomorphic graph up to edges ordering form graphqnet_list.

    :nqubit: int, the number of qubits
    :graphqnet_list:list(GraphQNet), list of the GraphQNet objects.
    :net_edges:list(int), list of all possible edges in network (integer) format
    """

    #for each net, spawn new unique net with one more edge
    unique_net_all = []
    for i,gqn in enumerate(graphqnet_list):
        unique_net = []
        for net in net_edges:
            gqn_cand = GraphQNet(nqubit, (*gqn.netgates, net))
            if not gqn_cand.more_three_con_edges():
                if not gqn_cand.is_isomorphic_uptolist(unique_net):
                    unique_net.append(gqn_cand)
        unique_net_all.extend(unique_net)

    return unique_net_all



def graphqnet_noniso(nqubit, net_depth, outdir=False, start_gqns=False, draw_graphs=False, conjugation_by_swap=True, time_reversal=False):
    """ List uninque non-isomorphic graph by iterating it

    :nqubit: int, the number of qubits
    :net_depth: int, the depth target
    :outdir: str='out' or boolean, the directory to store outputs. The output files
             will have format 'net-[nqubit]q-[nedges]e.json'
             that will store the unique 2-bit gates network.
    :start_gqns:list(GraphQNet), the list of unique GraphQNet object as starting point of iteration
    :draw_graphs:boolean, if draw all the resulting graphs. It will drawn inside the outdir folder
    :conjugation_by_swap: boolean=True, consider elimination by swap conjugation
    :time_reversal: boolean=True, consider elimination by time reversal
    """
    # setup directories, files, and initial variables
    outdir = outdir if outdir else 'out'
    run(['mkdir','-p',outdir])

    if start_gqns :
        nedge, gqn_list = start_gqns[0].depth, start_gqns
    else :
        nedge, gqn_list = 1, [GraphQNet(nqubit, (bitop.pos_ones_toint(0,1),))]

    net_edges = [bitop.pos_ones_toint(*e) for e in combinations(range(nqubit),2)]

    # iteration part
    while nedge < net_depth:
        start_time = time()
        gqn_list = iterate_graphqnet_noniso(nqubit, gqn_list, net_edges)
        nedge += 1

        #eliminate the conjugation by swaps
        if conjugation_by_swap:
            for i,gqn in enumerate(gqn_list):
                equiv_gqn = gqn.conjugation_by_swap()
                to_elim = False
                for gqn2 in gqn_list[i+1:]:
                    if gqn2.is_isomorphic_uptolist(equiv_gqn):
                        to_elim = True
                        break
                if to_elim :
                    gqn_list[i]=False
            gqn_list = [gq for gq in gqn_list if gq]

        #eliminate the time reversal
        if  time_reversal:
            for i,gqn in enumerate(gqn_list):
                equiv_gqn = gqn.equivnet_time_reversal()
                to_elim = False
                for gqn2 in gqn_list[i+1:]:
                    if gqn2.is_isomorphic_to(equiv_gqn):
                        to_elim = True
                        break
                if to_elim :
                    gqn_list[i]=False
            gqn_list = [gq for gq in gqn_list if gq]

        # storing results
        res_path = '%s/net-%iQ-%iE.json'%(outdir,nqubit,nedge)
        res = {'nqubit':nqubit,
               'time': time()-start_time,
               'conjugation_by_swap': conjugation_by_swap,
               'time_reversal': time_reversal,
               'networks':[gqn.netgates for gqn in gqn_list]
               }
        with open(res_path, 'w+') as outf :
            json.dump(res, outf)

        if draw_graphs :
            draw_path = '%s/net-%iQ-%iE.png'%(outdir,nqubit,nedge)
            if len(gqn_list) > 0 :
                GraphQNet.draw_netgraphs_list(res['networks'], nqubit, outfile=draw_path)
            else : print("empty result, no image is produced")

    return gqn_list





def unique2net(nqubit, net_depth, startfile=False, draw_graphs=True, dirpath='out', ncpu=False, conjugation_by_swap=True, time_reversal=False):
    """
    Get a list of 2-bit gates networks. The unique gates are iterated by the following steps:
        1) iterate the non-isomorphic graph up to gate ordering
        2) apply more criteria of DiVincenzo and Smolin: time reversal

    param
        :nqubit: int, the number of qubits
        :net_depth: int, the depth of the gate-networks
        :startfile: str=False, the file that stores unique networks. The file must
                    be in format of dictionary {'nqubits':int, 'networks':[(int,int),..]}.If
                    it is present, the iteration will be started from there.
        :draw_graphs: boolean=True, to draw the produced graphs
        :dirpath: str=out, directory path to store outputs
        :ncpu: int=cpu_count, the number of cpu in parallelization
        'conjugation_by_swap'
        :time_reversal: boolean=False, include time reversal criteria

    return
        [(int,int,..),(...),...] a list of 2-bit networks gates, with the LSB
        convention. For example:
          network (3,5) is
            q0 -o-o--
                | |
            q1 -o-|--
                  |
            q2 ---o--
        where (q0,q1)=3, (q0,q2)=5
    """
    ncpu = ncpu if ncpu else cpu_count()
    start_gqns = False
    if startfile:
        with open(startfile) as iff:
            gqn_list = json.load(iff)['networks']
            start_depth = len(gqn_list[0])
            start_gqns = [ GraphQNet(start_depth, x) for x in gqn_list]
            if start_depth >= net_depth :
                print("nothing to do here")
                return

    start=time()

    unique_net = graphqnet_noniso(nqubit, net_depth, outdir=dirpath,
                                  start_gqns=start_gqns, draw_graphs=draw_graphs,
                                  conjugation_by_swap=conjugation_by_swap, time_reversal=time_reversal)

    print('%i unique networks is calculated in %f seconds'%(len(unique_net),time()-start))

    return unique_net



