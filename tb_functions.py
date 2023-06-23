import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

def tb0_model(N,t,e):     
	''' Tight Binding (spinless) model for periodic boundary condition
	Input : Number of lattice sites, hopping constant (int/float), onsite energies (list)
	Output : Tight binding Hamiltonian, eigenvalues and eigenvectors of the matrix ''' 
    H = np.zeros((N, N))
    np.fill_diagonal(H,e)                 #filling on site energies in the diagonal of H

    pairs = [[i,(i+1)%N] for i in range(N)]
    for pair in pairs :                    #filling hopping terms in H
        H[pair[0]][pair[1]] = t
        H[pair[1]][pair[0]] = t
    eigval,eigvec = np.linalg.eigh(H)
    new_vec = list(zip(*eigvec))
    
    return H, eigval,new_vec

def tb0_state_occ(states) :
	'''Describes occupation of eigenstates over the lattice sites 
	Input : list of eigenstates
	Output : plots the occuptation number graph and 
	returns list of occupation number over sites for all thes states (list of lists)'''
    num_states = len(states)
    dim = len(states[0])
    states_occ = []
    for index in range(num_states) : 
        mod_sq = np.zeros(dim)
        for i in range(dim) : 
            mod_sq[i] = np.abs(list_state[index][i])**2
        states_occ.append(mod_sq)

        plt.plot(range(dim), mod_sq, label='vec {:.2f}'.format(index))
    plt.xlabel('sites')
    plt.ylabel('<n>')
    plt.title("Occupation of states wrt sites")
    plt.legend()

    return states_occ