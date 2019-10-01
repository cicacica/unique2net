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
unique2net(
    nqubit,
    net_depth,
    startfile=False,
    draw_graphs=True,
    dirpath='out',
    ncpu=False,
    conjugation_by_swap=True,
    time_reversal=False,
)
```
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
    ```
      network (3,5) is
        q0 -o-o--
            | |
        q1 -o-|--
              |
        q2 ---o--
    where (q0,q1)=3, (q0,q2)=5
    ```
