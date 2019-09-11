# unique2net

A module to list all unique 2-bit gate networks using the criteria of DiVincenzo and Smolin in 
reference [Results on two-bit gate design for quantum computers](https://arxiv.org/abs/cond-mat/9409111). 

The unique gates are iterated by the following steps:

- Step 1: 
    Iterate non-isomorphic graphs by adding edges one-by-one to a fixed number of nodes. 
    It has only a serial implementation, but it can start with a list of non-isomorphic graphs that 
    have fewer edges.
- Step 2:
    Assign the edges ordering to the non-isomorphic graphs from Step 1. 
    At this step, the relabelling qubit has already implemented.
- Step 3:
    Apply more criteria of Divincenzo and Smolin: time reversal and conjugation by swapping. 


### Get a list of 2-bit gates network
```sh
from unique2net import unique2net

L = unique2net(number_of_qubits, depth_of_networks)
```

#### With default setting
```sh
unique2net(nqubit, net_depth, path_json=False, draw_graphs=True, save_edges=True,
            dirpath=[out-number_of_qubits]qubit-depth[depth_of_networks],
            ds_bit_permutation=False, ds_conjugation_by_swap=True,
            ds_time_reversal= False, ncpu=cpu_count())
```

### Draw the graphs based on equivalence classes defined in DS paper
```sh
test_draw_by_equivalents(number_of_qubits,depth_of_networks,ncpu=4,images_per_row=4,dirpath='out')
```
the result for 4 qubits with depth 4 is shown in the folder example
