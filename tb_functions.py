import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh

def tb0_model(N,t,e): 
    ''' One-orbital Tight Binding (spinless) model for periodic boundary condition, 1 electron problem
    Input : Number of lattice sites (int), hopping constant (int/float), onsite energies (list)
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

def spinless_state_occupation(states) :
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
            mod_sq[i] = np.abs(states[index][i])**2
        states_occ.append(mod_sq)

        plt.plot(range(dim), mod_sq, label='vec {:.2f}'.format(index))
    plt.xlabel('sites')
    plt.ylabel('<n>')
    plt.title("Occupation of states wrt sites")
    plt.legend()

    return states_occ

def density_of_states(eigval,n_centers) : 
    '''
    Input : eigenvalues of the model, number of points where density is to be found
    Output : density of states, plots density of states graph 
    '''
    gamma = (eigval[-1]-eigval[0])/n_centers
    centers = np.linspace(eigval[0]-n_centers/10, eigval[-1]+n_centers/10, num=n_centers)
    
    density = np.zeros(n_centers)
    for i, center in enumerate(centers):
        if center < eigval[0] or center > eigval[-1]: 
            density[i] = 0 
        else : 
            lorentz = np.sum(1 / np.pi * (gamma / 2) / ((center - eigval)**2 + (gamma / 2)**2))
            density[i] = lorentz
    norm_density = [float(i)/sum(density) for i in density]
    plt.plot(centers,norm_density)
    plt.xlabel('Energy of system')
    plt.ylabel('Density of states')
    plt.title("Density of states")
    return density