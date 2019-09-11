# unique2net

A module to list all unique 2-bit gate networks using the criteria of DiVincenzo and Smolin in reference arXiv:cond-mat/9409111. 

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


### Usage
```sh
from unique2net import unique2net

L = unique2net(number_of_qubits, depth_of_networks)
```

### With default setting
```sh
unique2net(nqubit, net_depth, path_json=False, draw_graphs=True, save_edges=True,
            dirpath=[out-number_of_qubits]qubit-depth[depth_of_networks],
            ds_bit_permutation=False, ds_conjugation_by_swap=True,
            ds_time_reversal= False, ncpu=cpu_count())
```
