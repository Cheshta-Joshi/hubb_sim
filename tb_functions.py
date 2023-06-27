import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh
import math 
from itertools import combinations
from tabulate import tabulate

def spinless_sub_dim(N,r) : 
    ''' Subspace dimension for spinless models 
    Input : Number of lattice sites, number of electrons
    Output : dimension of subspace'''
    return math.comb(N,r)

def spinless_fock_dim(N) : 
    ''''Fock dimension for spinless models
    Input : Number of lattice sites
    Output : Dimension of fock space'''
    dim_dict = {} 
    for r in range(N+1) : 
        dim_dict[r] = spinless_sub_dim(N,r)
    return sum(dim_dict.values()), dim_dict

def spinless_basis(N,r) : 
    '''Basis set for spinless models
    Input : Number of lattice sites, number of electrons
    Output : Set of all basis states'''
    basis_set = []
    lattice = list(range(N))
    places = list(combinations(lattice, r))
    for combination in places : 
        basis = [False] *N
        for index in combination : 
            basis[index] = True 
        basis_set.append(basis)
    return basis_set

def tb0_model(N,r,t,e): 
    ''' Generalised Tight Binding (spinless) model for periodic boundary condition
    Input : Number of lattice sites (int), number of electrons(int), hopping constant (int/float), onsite energies (list)
    Output : Tight binding Hamiltonian, eigenvalues and eigenvectors of the matrix ''' 
    dim = spinless_sub_dim(N,r)
    H = np.zeros((dim, dim))
    basis_set = spinless_basis(N,r)
    n_diag = np.zeros(dim)
    for i,basis in enumerate(basis_set) : 
        n_diag[i] = np.count_nonzero(basis)
    np.fill_diagonal(H,n_diag*e)
    for basis_index,basis in enumerate(basis_set) : 
        for site in range(len(basis)) : 
            if basis[site] == False and basis[(site+1)%N] == True : 
                new_state = basis.copy()
                new_state[site] = True
                new_state[(site+1)%N] = False 
                for i in range(len(basis_set)) : 
                    if basis_set[i] == new_state: 
                        f_index = i
                H[f_index][basis_index] +=t

            if basis[site] == True and basis[(site+1)%N] == False : 
                new_state = basis.copy()
                new_state[site] = False
                new_state[(site+1)%N] = True 
                for i in range(len(basis_set)) : 
                    if basis_set[i] == new_state : 
                        f_index = i

                H[f_index][basis_index] +=t  
    eigval,eigvec = np.linalg.eigh(H)
    new_vec = list(zip(*eigvec))
                
    return H,eigval,new_vec

def spinless_state_vec(basis) : 
    ''' Vector Representation of a basis state in spinless models
    Input : basis (boolean list/array)
    Output : vector form of the state '''
    N = len(basis) 
    r = np.count_nonzero(basis)
    dim = spinless_sub_dim(N,r)
    basis_set = spinless_basis(N,r)
    i =-1
    index=0
    for state in basis_set : 
        i+=1
        vec = np.zeros(dim)
        if state == basis : 
            index = i

        vec[index] = 1
    return vec

def state_occupation(states) :
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
    centers = np.linspace(eigval[0], eigval[-1], num=n_centers)
    
    density = np.zeros(n_centers)
    for i, center in enumerate(centers):
        lorentz = np.sum(1 / np.pi * (gamma / 2) / ((center - eigval)**2 + (gamma / 2)**2))
        density[i] = lorentz
        
    plt.plot(centers,density)
    plt.xlabel('Energy of system')
    plt.ylabel('Density of states')
    plt.title("Density of states")
    return density  