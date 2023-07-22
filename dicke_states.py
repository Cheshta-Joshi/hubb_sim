import numpy as np 
from qiskit import * 
import math
from qiskit.circuit.library import RYGate
from qiskit.circuit import Parameter

#Creates circuit for nDk states 
def scs(n,l) : 
    circ = QuantumCircuit(n)
    circ.cx(-2,-1)
    circ.cry(2*np.arccos(np.sqrt(1/n)),-1,-2)
    circ.cx(-2,-1)
    
    for i in range(l-1) : 
        circ.cx(-3-i,-1)
        ccry = RYGate(2*np.arccos(np.sqrt((2+i)/n))).control(2,label=None)
        circ.append(ccry,[-1,-2-i,-3-i])
        circ.cx(-3-i,-1)
    return circ 

def dicke_circ(n,k) : 
    pairs = []
    for a in range(n,k,-1) : 
        pairs.append([a,k])
    for a in range(k,1,-1) : 
        pairs.append([a,a-1])

    dk_circ = QuantumCircuit(n)
    dk_circ.x(range(-1,-k-1,-1))
    for pair in pairs : 
        dk_circ.append(scs(pair[0],pair[1]),range(pair[0]))
    return dk_circ

#Use the following code to get the statevector of a circuit 
def sv(circ) : 
	backend = Aer.get_backend('statevector_simulator')
	job = execute(circ, backend)
	result = job.result()
	statevector = result.get_statevector()
	return statevector

#Creates parametrized circuit for nDk 
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
    for q in range(k) : 
        dk_circ.x(-(k+1))
    var=0
    for pair in pairs : 
        new_circ,new_var = scs_param(pair[0],pair[1],var,param)
        var = new_var
        dk_circ.append(new_circ, range(pair[0]))
    return dk_circ

