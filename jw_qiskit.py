import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, Pauli, SparsePauliOp
from qiskit_nature.second_q.hamiltonians.lattices import LineLattice, BoundaryCondition

N = 3
e = [1]*N
t = 1
H = np.zeros((2**N,2**N))
X = [[0,1],
    [1,0]]
Y = [[0,-1j],
    [1j,0]]
Plus = np.array([[0,1],
       [0,0]])
Minus = np.array([[0,0],
        [1,0]])

H_T = np.zeros((2**N,2**N))
for j in range(N) : 
    mat_t = Minus @ Plus 
    t1 = np.kron(np.eye(2**j),mat_t)
    t_f = np.kron(t1,np.eye(2**(N-j-1)))
    H_T += e[j]*t_f

pairs = []
for i in range(N) : 
    pairs.append([i,(i+1)%N])  
    
H_hopp = np.zeros((2**N,2**N))
for pair in pairs : 
    i = pair[0]
    j = pair[1]
    mat_list = [np.eye(2)] * N
    mat_list[i] = Minus
    mat_list[j] = Plus
    #print(mat_list)
    A = mat_list[0]
    for k in range(N-1) : 
        B = mat_list[k+1]
        tensor = np.kron(A,B)
        A = tensor 
        
    H_hopp += t*tensor
if N != 2 :   
    for pair in pairs : 
        i = pair[0]
        j = pair[1]
        mat_list = [np.eye(2)] * N
        mat_list[i] = Plus
        mat_list[j] = Minus
        #print(mat_list)
        A = mat_list[0]
        for k in range(N-1) : 
            B = mat_list[k+1]
            tensor = np.kron(A,B)
            A = tensor 

        H_hopp += t*tensor
H=H_hopp+H_T     
print(H)
    
#Sparse Pauli Operator in qiskit for the JW hamiltonian

op = SparsePauliOp.from_sparse_list([("I", [0], 1)], num_qubits=N)
for i in range(N) : 
    op += SparsePauliOp.from_sparse_list([("I", [i], 0.5),("Z", [i], -0.5),("XX", [i,(i+1)%N], 0.5),("YY", [i,(i+1)%N], -0.5)], num_qubits=N)

#Using qiskit's lattice model to get tight binding hamiltonian for verification
num_nodes = 4
boundary_condition = BoundaryCondition.PERIODIC
t = 5
e =1
line_lattice = LineLattice(
    num_nodes=num_nodes,
    edge_parameter=t,
    onsite_parameter=e,
    boundary_condition=boundary_condition,
)
set(line_lattice.graph.weighted_edge_list())

mat = line_lattice.to_adjacency_matrix(weighted=True)
print(tabulate(mat,tablefmt='plain'))