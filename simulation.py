#!/usr/bin/env python3

import numpy as np
from scipy.optimize import minimize




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
            if U(i,j) != 0.0 :
                objective_val += np.abs(U(i,j)-S_tot(i,j))

    return objective_val


        

def objective_function_big():
    pass




def main():
    
    # test 1
    # find the circuit for --H--
    #                      --H--
    #
    U1 = np.kron(Hadamard(),Hadamard())
    
    
    # initialize the guess X={x1,x2}, [0,2pi) as a list
    X = np.random.uniform(low=0.0, high=2*np.pi, size=2)

    result = minimize(objective_function, X, args=(U1), method="BFGS") 
    
    #should have result.val, result.x, and other

    #Check the result

    print(template_J_CPHASE(result.x))
    print(template_J_CPHASE(result['x']))

    




if __name__ == "__main__" :
    main()

