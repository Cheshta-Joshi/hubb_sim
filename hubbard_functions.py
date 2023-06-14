'''These functions are for an N site open lattice model with onsite energy, hopping and nearest interaction terms, '''
import numpy as np 
import math 
from itertools import combinations
import matplotlib.pyplot as plt
from scipy.linalg import eigh, block_diag
from tabulate import tabulate
import time

def subspace_dim(N,r) : 
    '''
    input : number of lattice sites (N), number of electrons (r) 
    output : dimension of this subspace 
    '''
    return math.comb(N,r)

def fock_dim_dict(N) : 
    '''
    input : N => Number of lattice points
    output : Dimension of fock space (dim_fock), dictionary of all subspaces with dimension { (r) : dim_sub } 
    '''
    dim_dict = {}
    for r in range(N+1) : 
        dim_dict[r] = subspace_dim(N,r)
    return sum(dim_dict.values()), dim_dict

def basis_set(N,r) : 
    '''
    input = number of lattice sites (N), number of electrons (r) 
    output = list of basis [0110] example for N=4, r=2
    '''
    basis_dim = subspace_dim(N,r)
    basis_set = np.zeros((basis_dim,basis_dim), dtype=bool)
    choices = list(combinations(list(range(N)), r))
    for i in range(len(choices)) : 
        basis = np.zeros(basis_dim, dtype=bool)
        for index in choices[i] : 
            basis[index] = True
        basis_set[i] = basis
        
def H_subspace(N,r,e,t,U) :
    '''
    input = N =Number of lattice sites, r = number of electrons, 
    e=epsilon/onsite energy, t= hopping constant, U= interaction coefficient
    uses function subspace_dim, basis_set
    output = Hamiltonian of subspace
    '''
    sub_dim = subspace_dim(N,r)
    H_sub = np.zeros((sub_dim,sub_dim))
    basis = basis_set(N,r)
    #H_D
    np.fill_diagonal(H_sub, e*r)

    #H_T
    H = np.zeros((sub_dim,sub_dim))
    index = -1
    for state in basis : 
        index += 1
        final_list = []
        final_index = []
        for i in range(len(state)-1) : 
            if state[i] == False and state[i+1] == True  :
                new_state = state.copy()
                new_state[i] = True
                new_state[i+1] = False 
                final_list.append(new_state)
                final_index.append(np.where((basis == new_state).all(axis=1))[0][0])
            if state[i] == True and state[i+1] == False :
                new_state = state.copy()
                new_state[i] = False
                new_state[i+1] = True 
                final_list.append(new_state)
                final_index.append(np.where((basis == new_state).all(axis=1))[0][0])
        #print(final_index)
        if final_index != [] :
            #print("yes")
            for f in final_index : 
                #print(f)
                #print(index)
                H_sub[f][index] += 1*t
    #H_U
    U_diag = []
    for state in basis :
        count = 0
        for i in range(N-1) : 
            if state[i] == True and state[i+1] == True : 
                count +=1
        U_diag.append(U*count)
    H_sub += np.diag(U_diag)
    return H_sub

def full_hamiltonian(N,e,t,U): 
    ''' 
    input = N,e,t,u; uses function H_subspace
    output = Hamiltonian of full fock space
    '''
    H_sub_list = []
    for r in range(N+1) : 
        H_sub = H_subspace(N,r,e,t,U)
        print("H_sub:\n",H_sub)
        e_sub, v_sub = np.linalg.eigh(H_sub)
#         print("Eigenvalues of H_sub : \n",e_sub)
#         print("Eigenvectors of H_sub : \n",v_sub)
        print("---------")
        #print(e_sub)
        H_sub_list.append(H_sub)
    H_fock = block_diag(*H_sub_list) 
    return H_fock

''' Below functions are for an N site open lattice model, where electrons can be in up or down spin 
and the hamiltonian includes the hopping and onsite energy terms ''' 

def s2_sub_dim(N,r_up, r_down) : 
    '''
    input : number of lattice sites (N), number of spin up (r_up) and spin down electrons (r_down)
    output : dimension of this subspace 
    '''
    dim = math.comb(N,r_up) * math.comb(N,r_down)
    return dim

def s2_fock_dim_dict(N) : 
    '''
    input : N => Number of lattice points
    output : Dimension of fock space (dim_fock), dictionary of all subspaces with dimension { (r_up,r_down) : dim_sub } 
    '''
    dim_dict = {}
    pairs = []  # different combinations of number of up and down spins for fock space
    for n_up in range(N+1) : 
        for n_down in range(N+1): 
            pairs.append((n_up,n_down))
    dim_fock = 0 
    i = -1
    for r_up, r_down in pairs : 
        i+=1
        sub_dim = s2_sub_dim(N,r_up,r_down)
        dim_dict[pairs[i]] = sub_dim
        dim_fock += sub_dim
    return dim_fock, dim_dict
    
def s2_basis_set(N,r_up,r_down) : 
    '''
    input = number of lattice sites (N), number of spin up (r_up) and spin down electrons (r_down)
    output = list of basis [0110][0001] examp-le for N=4, r_up = 2, r_down = 1
    '''
    sub_arrays = []
    lattice = list(range(N))
    up_choices = list(combinations(lattice, r_up))
    for up_set in up_choices : 
        up = [False]*N
        for index in up_set : 
            up[index] = True
        down_choices = list(combinations(lattice, r_down))
        i = -1
        for down_set in down_choices : 
            down = [False]*N
            for index in down_set : 
                down[index] = True
            sub_arrays.append(np.array([up,down]))
            i+=1
    basis = np.array(sub_arrays)
    return basis

def s2_full_hamiltonian(N,t,U) :
    basis = s2_basis_set(N,r_up,r_down)
    sub_dim = len(basis)
    hubbard = np.zeros((sub_dim,sub_dim))
    true_counts = []
    for state in basis : 
        el_and = [a and b for a, b in zip(state[0], state[1])]
        true_counts.append(el_and.count(True))  
    diag = [U*t for t in true_counts ]
    np.fill_diagonal(hubbard,diag)
    index = -1
    for state in basis : 
        index+=1
        final_list = []
        final_index = []
        for i in range(N-1) :
            for spin in range(2) : 
                if state[spin][i] == False and state[spin][i+1] == True : 
                    new_state = state.copy()
                    new_state[spin][i] = True 
                    new_state[spin][i+1] = False 
                    final_list.append(new_state)
                    final_index.append(np.where(np.all(basis == new_state, axis=(1, 2)))[0][0])
                if state[spin][i] == True and state[spin][i+1] == False : 
                    new_state = state.copy()
                    new_state[spin][i] = False  
                    new_state[spin][i+1] = True 
                    final_list.append(new_state)
                    final_index.append(np.where(np.all(basis == new_state, axis=(1, 2)))[0][0])
        if final_index != [] :
                #print("yes")
                for f in final_index : 
                    hubbard[f][index] += 1*t
    return(hubbard)