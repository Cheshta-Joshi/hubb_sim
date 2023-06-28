'''These functions are for an N site open lattice model with onsite energy, hopping and nearest interaction terms, '''
import numpy as np 
import math 
from itertools import combinations
import matplotlib.pyplot as plt
from scipy.linalg import eigh, block_diag
from tabulate import tabulate

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
        
def hubb0_model(N,r,t,e,U): 
    ''' Generalised Tight Binding (spinless) model for periodic boundary condition
    Input : Number of lattice sites (int), number of electrons(int), hopping constant (int/float), onsite energies (list), interaction term U (int)
    Output : Tight binding Hamiltonian, eigenvalues and eigenvectors of the matrix ''' 
    dim = spinless_sub_dim(N,r)
    #Special Cases
    if r==0 : 
        H = np.zeros(1)
        eigval = H[0]
        new_vec = H
    elif dim == 1 and r!= 0 : 
        H = np.array(e[0]+U)
        eigval = H
        new_vec = H
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
        for i in range(dim) :
            count = 0
            for j in range(N) : 
                if basis_set[i][j] == 1 and basis_set[i][(j+1)%N] == 1: 
                    count += 1 
            H[i][i] += count * U
        eigval,eigvec = np.linalg.eigh(H)
        new_vec = list(zip(*eigvec))                   
    return H,eigval,new_vec

def hubb0_full(N,e,t,U): 
    '''Full Block Diagonal Hamiltonian for some N length closed lattice '''
    H_sub_list = []
    for r in range(N+1) : 
        H_sub,e_sub,v_sub = hubb0_model(N,r,t,e,U)
        H_sub_list.append(H_sub)
    H = block_diag(*H_sub_list) 
    return H

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