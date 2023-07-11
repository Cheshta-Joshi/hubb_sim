import numpy as np
from qiskit import *
from qiskit.quantum_info import Operator, Pauli, SparsePauliOp
from tabulate import tabulate
from itertools import combinations

def spinless_basis(N,r) : 
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

def ferm_JW(mat) : 
    '''Input : SparsePauliOp JW matrix, uses function spinless_states_index
    Output : Returns JW matrix in order of fermionic basis states'''
    class_index = spinless_states_index(N)
    JW_mat = mat.to_matrix()
    ferm_mat = np.zeros((2**N,2**N))
    for i in range(2**N) : 
        for j in range(2**N) : 
            ferm_mat[i][j] = JW_mat[class_index[i]][class_index[j]]
    return ferm_mat
