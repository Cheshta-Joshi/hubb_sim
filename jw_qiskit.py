import numpy as np
from qiskit import QuantumCircuit
from qiskit.quantum_info import Operator, Pauli, SparsePauliOp

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
    
