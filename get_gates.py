#!/usr/bin/env python3

from unique2net import unique2net, all_2g
import json
import sys



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
                            

                        

