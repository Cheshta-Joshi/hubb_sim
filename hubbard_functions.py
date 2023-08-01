'''These functions are for an N site open lattice model with onsite energy, hopping and nearest interaction terms, '''
import numpy as np 
import math 
from itertools import combinations
import matplotlib.pyplot as plt
from scipy.linalg import eigh, block_diag
from tabulate import tabulate

#Common functions
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

#Hubbard model functions      
def hubb0_model(N,r,e,t,U): 
    ''' Generalised Tight Binding (spinless) model for periodic boundary condition
    Input : Number of lattice sites (int), number of electrons(int), hopping constant (int/float), onsite energies (list), interaction term U (int)
    Output : Tight binding Hamiltonian, eigenvalues and eigenvectors of the matrix ''' 
    dim = spinless_sub_dim(N,r)
    #Special Cases
    if r==0 : 
        H = np.zeros(1)
        eigval = 0
        new_vec = [[1]]
    elif r==N : 
        H = [[sum(e)+N*U]]
        eigval = H[0]
        new_vec = [[1]]
    if N == 1 and r==1: 
        H = [e]
        eigval = H[0]
        new_vec = [[1]]
    else : 
        H = np.zeros((dim, dim))
        basis_set = spinless_basis(N,r)
        #H_diagonal, onsite energy terms
        n_diag = np.zeros(dim)
        for i in range(dim) : 
            for j in range(N) : 
                n_diag[i] += e[j]*basis_set[i][j]
        np.fill_diagonal(H,n_diag)
        #H_T Hopping terms 
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
        #H_U, interaction terms
        for basis_index,basis in enumerate(basis_set) : 
            for site in range(len(basis)) : 
                if basis[site] == True and basis[(site+1)%N] == True : 
                    H[basis_index][basis_index] +=U

        eigval,eigvec = np.linalg.eigh(H)
        new_vec = list(zip(*eigvec))                   
    return H,eigval,new_vec

def hubb0_full(N,e,t,U): 
    '''Full Block Diagonal Hamiltonian for some N length closed lattice '''
    H_sub_list = []
    for r in range(N+1) : 
        H_sub = hubb0_model(N,r,e,t,U)[0]
        H_sub_list.append(H_sub)
    H = block_diag(*H_sub_list) 
    return H

#JW functions

def hubb0_JW(N,e,t,U) : 
    strings = []
    opt = SparsePauliOp.from_sparse_list([("I", [0], 0)], num_qubits=N)  
    for k in range(N) : 
        a0='I'*N
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
        
        c_base='I'*N

        c0 = c_base
        
        c1_list = list(c_base)
        c1_list[k] = 'I'
        c1_list[(k+1)%N] = 'Z'
        c1 = ''.join(c1_list)
        
        c2_list = list(c_base)
        c2_list[k] = 'Z'
        c2_list[(k+1)%N] = 'I'
        c2 = ''.join(c2_list)
        
        c3_list = list(c_base)
        c3_list[k] = 'Z'
        c3_list[(k+1)%N] = 'Z'
        c3 = ''.join(c3_list)
        T=t
        if N==2 : 
            T=t/2
        
        opt += SparsePauliOp.from_list([(a0, 0.5*e[k]), (a1, -0.5*e[k]),
                                        (new_b0, 0.5*T),(new_b1, 0.5*T),
                                        (c0, U*0.25),(c1, -0.25*U),(c2, -0.25*U),(c3,U*0.25)])
    return opt

#Dicke states
def scs_param(n,l,var,param) : 

    circ = QuantumCircuit(n)
    circ.cx(-2,-1)
    circ.cry(param[var],-1,-2)
    var+=1
    circ.cx(-2,-1)

    for i in range(l-1) : 
        circ.cx(-3-i,-1)
        ccry = RYGate(param[var]).control(2,label=None)
        var+=1
        circ.append(ccry,[-1,-2-i,-3-i])
        circ.cx(-3-i,-1)
    return circ, var

def dicke_param(n,k) : 
    pairs = []
    for a in range(n,k,-1) : 
        pairs.append([a,k])
    for a in range(k,1,-1) : 
        pairs.append([a,a-1])

    num_angles = int(k*(n-k) + k*(k-1)/2)
    param = [Parameter(f"angle_{i+1}") for i in range(num_angles)]

    dk_circ = QuantumCircuit(n)
    for ind in range(k) : 
        dk_circ.x(-(1+ind))
    var=0
    for pair in pairs : 
        new_circ,new_var = scs_param(pair[0],pair[1],var,param)
        var = new_var
        dk_circ.append(new_circ, range(pair[0]))
    return dk_circ

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

#properties
def spinless_state_vec(basis) : 
    '''Vector represnetation of basis state
    Input : Boolean list/array of a basis
    Output : Vector representation of the basis'''
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
    
def spinless_state_occupation(states) :
    '''Describes occupation of eigenstates over the lattice sites 
    Input : list of eigenstates
    Output : plots the occuptation number graph and 
    returns list of occupation number over sites for all thes states (list of lists)'''
    num_states = len(states)
    dim = len(states[0])
    vec_site_occ = []
    basis_set = spinless_basis(N,r)
    for index,vec in enumerate(states) : 
        site_occ = np.zeros(N)
        for site in range(N) : 
            for i in range(dim) : 
                site_occ[site] += np.abs(vec[i])**2 *basis_set[i][site]
        plt.plot(range(N), site_occ, label='vec {:.2f}'.format(index))
        vec_site_occ.append(site_occ)
    plt.xlabel('sites')
    plt.ylabel('<n>')
    plt.title("Occupation of states wrt sites")
    plt.legend()
    return site_occ

def density_of_states(eigval,n_centers,zoom) :
    '''Energy density of states
    Input : eigenvalues (list), Number of points for density calculation, Zoom the density part of the plot
    Output : Plot of energy density of states]'''
    gamma = (eigval[-1]-eigval[0])/n_centers
    centers = np.linspace(eigval[0]-n_centers/zoom, eigval[-1]+n_centers/zoom, num=n_centers)
    
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