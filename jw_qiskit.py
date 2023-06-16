from qiskit import * 
import numpy as np
from qiskit.opflow import (I, X,Y,Z,StateFn, Zero, One, Plus, Minus, H, CX,
                           DictStateFn, VectorStateFn, CircuitStateFn, OperatorStateFn)

N = 4 
e = 1
t = 2
U = 3

initial_state = One^N
s_plus = (X - 1j*Y)/2
s_minus = (X + 1j*Y)/2

r = 1
op = e*((s_plus @ s_minus) ^ I) + U*((s_plus @ s_minus)^(s_plus @ s_plus)) + t *((s_plus ^ s_minus)+ (s_minus ^s_plus))

x = r
y = N-r-2
full_op = (I^x) ^ op ^ (I^y) 

print(full_op.eval(initial_state))



