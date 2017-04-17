import numpy as np

def int_to_bits(n, num_bits):
    a = np.ndarray(num_bits, np.bool)
    for i in range(num_bits):
        a[i] = (n >> i) & 1
    return a



def expand_gate(num_qubits, gate, to_qubits):
    N = 2**num_qubits
    to_qubits_set = set(to_qubits)
    M = np.zeros(N, N)
    for qubit in range(num_qubits):
        if qubit in to_qubits_set:

    for 

    for same_index in range(num_qubits - len(to_qubits)):
        same_index_string = int_to_bits(same_index)


        
