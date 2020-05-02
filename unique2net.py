#!/usr/bin/env python3

__doc__=""" unique2net.py: list all unique gates criteria of
DiVincenzo and Smolin (cond-mat/9409111).

The unique gates are iterated by the following steps:
    1) iterate the non-isomorphic graph up to ordering(non-parallell)
    2) At this step, the relabelling qubit has already implemented.
    3) apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping


MAIN USAGE:

    from unique2net import unique2net

    unique2net(nqubit, network_length)

"""
__author__ = "Cica Gustiani"
__license__ = "GPL"
__version__ = "4.1.0"
__maintainer__ = "Cica Gustiani"
__email__ = "cicagustiani@gmail.com"



#standard libraries
from itertools import combinations
from time import time
from subprocess import run
from multiprocessing import cpu_count, Pool
import json
import os

#additional library
from GraphQNet import GraphQNet, bitop




def iterate_graphqnet_noniso(nqubit, graphqnet_list, net_edges):
    """
    Produce non-isomorphic graph up to edges ordering form graphqnet_list.
    It applies criteria: bit-permutation and cojugation by swap

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



def __helper_idx_conjugation_by_swap(args):
    return __idx_conjugation_by_swap(*args)


def __idx_conjugation_by_swap(gqn_list, idx):
    """
    return the index when it is equivalent by swap conjugation
    :gqn_list: the whole GraphQNet list
    :idx: the index to be checked
    """
    equiv_gqn = gqn_list[idx].conjugation_by_swap()
    for gqn2 in gqn_list[idx+1:]:
        if gqn2.is_isomorphic_uptolist(equiv_gqn):
            return idx
    return False


def __helper_idx_time_reversal(args):
    return __idx_time_reversal(*args)


def __idx_time_reversal(gqn_list, idx):
    """
    return the index when it is equivalent by swap conjugation
    :gqn_list: the whole GraphQNet list
    :idx: the index to be checked
    """
    equiv_gqn = gqn_list[idx].time_reversal()
    for gqn2 in gqn_list[idx+1:]:
        if gqn2.is_isomorphic_to(equiv_gqn):
            return idx
    return False




def graphqnet_noniso(nqubit, net_depth, outdir=False, start_gqns=False, draw_graphs=False, ncpu=False, conjugation_by_swap=True, time_reversal=False):
    """ List uninque non-isomorphic graph by iterating it

    :nqubit: int, the number of qubits
    :net_depth: int, the depth target
    :outdir: str='out' or boolean, the directory to store outputs. The output files
             will have format 'net-[nqubit]q-[nedges]e.json'
             that will store the unique 2-bit gates network.
    :start_gqns:list(GraphQNet), the list of unique GraphQNet object as starting point of iteration
    :draw_graphs:boolean, if draw all the resulting graphs. It will drawn inside the outdir folder
    :ncpu:int=cpu_count(),the cpu number for parallelization
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
    ncpu = ncpu if ncpu else cpu_count()
    # iteration part
    while nedge < net_depth:
        start_time = time()

        #check if such a result exists
        res_path = '%s/nonisonet-%iQ-%iE.json'%(outdir,nqubit,nedge+1)

        try :
            #load from previous calculation
            with open(res_path) as inff :
                res = json.load(inff)
            res['networks']=[tuple(ng) for ng in res['networks']]
            nedge += 1
            gqn_list = [GraphQNet(nqubit, ng) for ng in res['networks']]

            print(gqn_list, res_path)

        except FileNotFoundError:
            #do everything 

            gqn_list = iterate_graphqnet_noniso(nqubit, gqn_list, net_edges)
            nedge += 1

            #eliminate the conjugation by swaps
            if conjugation_by_swap:
                lenl = len(gqn_list)
                P = Pool(ncpu)
                to_elim = P.map(__helper_idx_conjugation_by_swap, zip([gqn_list]*lenl, range(lenl)))
                P.close()
                P.join()

                for i in filter(lambda x: x, to_elim) :
                    gqn_list[i] = False
                gqn_list = [gq for gq in gqn_list if gq]

            # storing results
            res = {'nqubit':nqubit,
                   'depth':nedge,
                   'time': time()-start_time,
                   'conjugation_by_swap': conjugation_by_swap,
                   'time_reversal': False,
                   'start_gate': start_gqns[0].depth  if start_gqns else 1,
                   'networks':[gqn.netgates for gqn in gqn_list]
                   }
            with open(res_path, 'w+') as outf :
                json.dump(res, outf)

        if draw_graphs :
            draw_path = '%s/nonisonet-%iQ-%iE.png'%(outdir,nqubit,nedge)
            if os.path.exists(draw_path):
                print("No figure generated, it has already done")
            else : 
                if len(gqn_list) > 0 :
                    GraphQNet.draw_netgraphs_list(res['networks'], nqubit, outfile=draw_path)
                else : print("empty result, no image is produced")


    #final result paths
    fin_path = '%s/net-%iQ-%iE.json'%(outdir,nqubit,net_depth)
    fin_draw_path = '%s/net-%iQ-%iE.png'%(outdir,nqubit,net_depth)

    #eliminate the time reversal
    if time_reversal:
        lenl = len(gqn_list)
        P = Pool(ncpu)
        to_elim = P.map(__helper_idx_time_reversal, zip([gqn_list]*lenl, range(lenl)))
        P.close()
        P.join()
        for i in filter(lambda x: x, to_elim) :
            gqn_list[i] = False
        gqn_list = [gq for gq in gqn_list if gq]
        res = {'nqubit':nqubit,
               'depth':nedge,
               'time': time()-start_time,
               'conjugation_by_swap': conjugation_by_swap,
               'time_reversal': time_reversal,
               'start_gate': start_gqns[0].depth  if start_gqns else 1,
               'networks':[gqn.netgates for gqn in gqn_list]
               }
        with open(fin_path, 'w+') as outf :
            json.dump(res, outf)

        if draw_graphs :
            GraphQNet.draw_netgraphs_list(res['networks'], nqubit, outfile=fin_draw_path)

    else : 
        part_path = '%s/nonisonet-%iQ-%iE.json'%(outdir,nqubit,net_depth)
        os.rename(part_path,fin_path)

        if draw_graphs :
            draw_path = '%s/nonisonet-%iQ-%iE.png'%(outdir,nqubit,net_depth)
            os.rename(draw_path, fin_draw_path)


    return gqn_list





def unique2net(nqubit, net_depth, startfile=False, draw_graphs=True, outpath='out', ncpu=False, conjugation_by_swap=True, time_reversal=True):
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
        :outpath: str=out, directory path to store outputs
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
            sdata = json.load(iff)
            sdepth = sdata['depth']
            gqn_list = sdata['networks']
            snqubit = sdata['nqubit']
            start_gqns = [GraphQNet(snqubit, x) for x in gqn_list]
            if sdepth >= net_depth :
                print("nothing to do here")
                return

    start=time()

    unique_net = graphqnet_noniso(nqubit, net_depth, outdir=outpath,
                                  start_gqns=start_gqns, draw_graphs=draw_graphs,
                                  ncpu = ncpu, conjugation_by_swap=conjugation_by_swap, time_reversal=time_reversal)

    print('%i unique networks is calculated in %f seconds'%(len(unique_net),time()-start))

    return unique_net



