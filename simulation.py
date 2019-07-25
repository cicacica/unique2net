#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize
import networkx as nx
import os




### collections of quantum gates

def Hadamard():
    """ the hadamard gate """
    return np.matrix([[1,1],[1,-1]])/np.sqrt(2)



def gate_J(x):
    """ J gate = [[1, 1], [e^{ix},-e^{ix}]] 
    
    param:
        :x: float, the parameter of J 
    """
    return np.matrix([[1,1],[np.exp(complex(0,x)), -np.exp(complex(0,x))]])




def template_J_CPHASE(X):
    """
    a big function, complicated

    param:
        :X: list(float) the parameters of gate_J


    return:
        The S_tot matrix
    """
    
    #temporary example, works only for --- J ---
    #                                  --- J ---

    return np.kron(gate_J(X[0]), gate_J(X[1]))




#### definitions of objective functions
## objective function always return an objective value, a float!

def objective_function(U, S_tot):
    """ the small f

    param:
        :U: np.array, the desired unitary
        :S_tot: np.array, the guessed matrix from a network gates

    return: 
        float, objective value
    """
    objective_val = 0.0 

    dim = len(U)

    for i in range(dim):
        for j in range(dim):
            if U[i,j] != 0.0 :
                objective_val += np.abs(U[i,j]-S_tot[i,j])

    return objective_val


        

def objective_function_big():
    pass




def graph_from_flow_paper():
    """
    The graph from "Determinism in the one-way model"
    provide the positions yourself
    """
    G = nx.Graph()
    G.add_node(1, pos=(0,0.5))
    G.add_node(2, pos=(0,1.5))
    G.add_node(3, pos=(1,1.5))
    G.add_node(4, pos=(2,2  ))
    G.add_node(5, pos=(1,0.5))
    G.add_node(6, pos=(2,1))
    G.add_node(7, pos=(2,0))
    G.add_edges_from([(2,3),(2,5),(3,4),(1,5),(3,4),(3,6),(5,6),(5,7)])    

    return G


def __draw_edge(source, sink, pos_source, pos_sink, arrow=False):
    """
    Draw an edge that may be a directed or non-directed edge

    :source: str, the node name of source
    :sink: str, the node name of sink
    :pos_source: (int,int), the coordinate of sink vertex
    :pos_sink: (int,int), the coordinate of sink vertex
    :arrow: boolean, if it's true, an arrow will be drawn from source to
            sink
    """
    dx = abs(pos_source[0]-pos_sink[0])
    dy = abs(pos_source[1]-pos_sink[1])

    ops = '' #optional decorations
    if arrow : ops += 'postaction={decorate}'

    draw = '\\draw[%s]'%ops
    if dx==0 and dy>1 :
        draw+= '(%s) to[out=-50,in=50] (%s);'%(source,sink)
    else :
        draw+= '(%s) -- (%s);'%(source,sink)

    return draw











def graph_to_texpdf(G, **kwargs):
    """
    Convert graph to latex

    G -- graph, the plotted graph

    **kwargs :
        :pos: boolean(def=True), if true the node must present attribute 'pos'
        :total_order: boolean(def=True), if true the node must present attribute 'sorting_key'
        :outpath: str(def=plot.tex) the path of latex output file
    """
    opt = {'pos': True, 'outpath': 'plot.tex', 'total_order': False,
            'V_flow': False}

    for key in opt :
            if key in kwargs : opt[key]=kwargs[key]

    tex = ['\\documentclass{standalone}']
    tex.append('\\usepackage{tikz,graphicx}')
    tex.append('\\usetikzlibrary{decorations.markings}')
    tex.append('\\newcommand{\\midarrow}{\\tikz \\draw[-triangle 90] (0,0) -- +(.1,0);}')
    tex.append('\\begin{document}')
    tex.append('\\begin{tikzpicture}')
    tex.append('[font=\\small,scale=.7,auto,every node/.style={circle,inner sep=.5,minimum size=7.8pt,draw=black,text width=7.9pt,align=center}]')

    if not opt['pos']:     #assign random positions if they are not assigned
         assign_random_positions(G)

    node_names = dict([(n,'n%s'%(str(n))) for n in G.nodes])  #names for nodes

    for n in G.nodes :
        total_order =  str(G.nodes[n]['sorting_key'])  if opt['total_order'] else str(n)
        vflow = G.nodes[n]['flow_class'] if opt['V_flow'] else ''

        node_number = str(n)
        tex.append('\\node[%s] (%s) at (%i,%i){%s};'%('', node_names[n], G.nodes[n]['pos'][0], G.nodes[n]['pos'][1],node_number))

    for n1, n2 in G.edges:
        tex.append(__draw_edge(node_names[n1],node_names[n2],G.nodes[n1]['pos'],G.nodes[n2]['pos']))

    tex.append('\\end{tikzpicture}')
    tex.append('\\end{document}')

    with open(opt['outpath'],'w') as outf :
        outf.write('\n'.join(tex))

    os.system("pdflatex --file-line-error  --shell-escape --synctex=1 %s"%opt['outpath'])



def main():
    
    # test 1
    # find the circuit for --H--
    #                      --H--
    #
#   U1 = np.kron(Hadamard(),Hadamard())
#   
#   print(len(U1))
#   
#   # initialize the guess X={x1,x2}, [0,2pi) as a list
#   X = np.random.uniform(low=0.0, high=2*np.pi, size=2)

#   result = minimize(objective_function, X, args=(U1), method="BFGS") 
#   
#   #should have result.val, result.x, and other

#   #Check the result

#   print(template_J_CPHASE(result.x))
#   print(template_J_CPHASE(result['x']))

    # test 2
    # draw a graph state 
    G = graph_from_flow_paper()
    graph_to_texpdf(G)

    




if __name__ == "__main__" :
    main()

