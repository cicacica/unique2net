#!/usr/bin/env python3

from unique2net import test_draw_by_equivalents
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

    nqubit, net_depth = int(sys.argv[1]),int(sys.argv[2])
    test_draw_by_equivalents(nqubit, net_depth, ncpu=16, images_per_row=8,
            dirpath='out-%i-%i'%(nqubit,net_depth)) 
    

