import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh, block_diag
from qiskit import *
from qiskit.quantum_info import Operator, Pauli, SparsePauliOp
from tabulate import tabulate
from itertools import combinations
import math

#common functions
def spinless_sub_dim(N,r) : 
    '''
    input : number of lattice sites (N), number of electrons (r) 
    output : dimension of this subspace 
    '''
    return math.comb(N,r)

def spinless_fock_dim(N) : 
    '''
    input : N => Number of lattice points
    output : Dimension of fock space (dim_fock), dictionary of all subspaces with dimension { (r) : dim_sub } 
    '''
    dim_dict = {} 
    for r in range(N+1) : 
        dim_dict[r] = spinless_sub_dim(N,r)
    return sum(dim_dict.values()), dim_dict

def spinless_basis(N,r) : 
    '''
    input = number of lattice sites (N), number of electrons (r) 
    output = list of basis [0110] example for N=4, r=2
    '''
    basis_set = []
    lattice = list(range(N))
    places = list(combinations(lattice, r))
    for combination in places : 
        basis = [False] *N
        for index in combination : 
            basis[index] = True 
        basis_set.append(basis)
    return basis_set

def spinless_states_index(N) : 
    '''Input : N, uses spinless_basis function
    Output : Returns binary of basis states in fermionic Hamiltonian ''' 
    index = []
    for r in range(N+1) : 
        b_set = spinless_basis(N,r)
        new_bset = []
        for b in b_set : 
            new_b = sum([2**(N-i-1)*b[i] for i in range(N)])
            index.append(new_b)
    return index


#tight binding model functions
def tb0_model(N,r,t,e): 
    ''' Generalised Tight Binding (spinless) model for periodic boundary condition
    Input : Number of lattice sites (int), number of electrons(int), hopping constant (int/float), onsite energies (list)
    Output : Tight binding Hamiltonian, eigenvalues and eigenvectors of the matrix ''' 
    dim = spinless_sub_dim(N,r)
    if r==0 : 
        H = np.zeros(1)
        eigval = H[0]
        new_vec = H
    elif r==N : 
        H = [[r]]
        eigval = H[0]
        new_vec = H
    else : 
        H = np.zeros((dim, dim))
        basis_set = spinless_basis(N,r)
        n_diag = np.zeros(dim)
        for i in range(dim) : 
            for j in range(N) : 
                n_diag[i] += e[j]*basis_set[i][j]

        np.fill_diagonal(H,n_diag)
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
                if N != 2 : 
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

def tb0_full(N,e,t): 
    '''Full Block Diagonal Hamiltonian for some N length closed lattice '''
    H_sub_list = []
    for r in range(N+1) : 
        H_sub = tb0_model(N,r,t,e)[0]
        H_sub_list.append(H_sub)
    H = block_diag(*H_sub_list) 
    return H

#JW functions
def tb0_JW(N,e,t) : 
    strings = []
    opt = SparsePauliOp.from_sparse_list([("I", [0], 0)], num_qubits=N)  
    for k in range(N) : 
        a0='I'*(N)
        a1 = 'I'*(k)+'Z' +'I'*(N-k-1)

        b0='I'*N
        b0_list = list(b0)
        b0_list[k] = 'X'
        b0_list[(k+1)%N] = 'X'
        new_b0 = ''.join(b0_list)

        b1='I'*N
        b1_list = list(b0)
        b1_list[k] = 'Y'
        b1_list[(k+1)%N] = 'Y'
        new_b1 = ''.join(b1_list)

        strings.append([a0,a1,new_b0,new_b1])
        opt += SparsePauliOp.from_list([(a0, 0.5*e[k]), (a1, -0.5*e[k]),(new_b0, 0.5*t),(new_b1, 0.5*t)])
    return opt  

def ferm_JW(JW_mat) : 
    '''Input : SparsePauliOp JW matrix, uses function spinless_states_index
    Output : Returns JW matrix in order of fermionic basis states'''
    class_index = spinless_states_index(N)
    #JW_mat = mat.to_matrix()
    ferm_mat = np.zeros((2**N,2**N))
    for i in range(2**N) : 
        for j in range(2**N) : 
            ferm_mat[i][j] = JW_mat[class_index[i]][class_index[j]]
    return ferm_mat


#Properties of tb model
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