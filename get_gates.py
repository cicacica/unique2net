#!/usr/bin/env python3

from unique2net import unique2net, change_format, all_2g
import json
import sys

if __name__ == "__main__" :
    """
        Usage : get_gates.py nqubit net_depth ncpu
    """

    nqubit, net_depth, ncpu = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3])

    g2=all_2g(nqubit)
    G = unique2net(nqubit, net_depth, ncpu)

    fname='rijul_nqubit-%i_depth-%i.txt'%(nqubit, net_depth)
    with open(fname,'w+') as outfile:
         json.dump(change_format(g2,nqubit,G), outfile, indent=None)
    
    fname='net_nqubit-%i_depth-%i.json'%(nqubit, net_depth)
    with open(fname,'w+') as outfile:
         json.dump(G, outfile, indent=None)
                            

                        

