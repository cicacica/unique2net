# unique2net

A program to list all unique 2-bit gate networks using the criteria of DiVincenzo and Smolin (DS-criteria) in 
reference [Results on two-bit gate design for quantum computers](https://arxiv.org/abs/cond-mat/9409111). 

The unique gates are obtained by the following steps:
- **Step 1**
    Iterate non-isomorphic graphs by adding one edge to a fixed number of nodes. 
    The isomorphism counts also the order of placing the edge.
    It has only a serial implementation, but it can start with previous results that have less edges.
    At this step, the first DS-criteria has already implemented e.g., bit relabelling.
    
- **Step 2**
    Apply other DS-criteria: conjugation by swap and time reversal --- if necessary.
    
- **Step 3**
    Go to **Step 1** until reaching the desired number of edges  


### Main usage: to a list of 2-bit gates network
```sh
from unique2net import unique2net

L = unique2net(number_of_qubits, depth_of_networks)
```
## Start iteration from existing result
Here, calculate network with 5 qubits, depth 5, from the result of 5 qubits, depth 4.
```sh
from unique2net import unique2net

L = unique2net(5, 5, startfile='out/net-5Q-4E.json')
```

#### The default setting and docstring

```sh
unique2net(
    nqubit,
    net_depth,
    startfile=False,
    draw_graphs=True,
    dirpath='out',
    ncpu=cpu_count(),
    conjugation_by_swap=True,
    time_reversal=False,
)
```
<pre>
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
   :list(GraphQNet): a list of GraphQNet objects, where the 2-bit gates network are unique 
                     according to the given criteria
   
     </pre>

## Running example, the results are present in folder 'out' 
```sh
from unique2net import unique2net

L = unique2net(5, 5, time_reversal=True)
```
