import numpy as np
from itertools import product
from scipy import sparse

def int_to_bits(n, num_bits):
    a = np.empty(num_bits, np.bool)
    for i in range(num_bits):
        a[i] = (n >> i) & 1
    return a

def int_to_bit_string(n, num_bits):
    s = bin(n)[2:]
    return ("0" * (num_bits - len(s)) + s)[::-1]

def bit_vector_to_int(a):
    n = 0
    for i, bit in enumerate(a):
        if bit:
            n |= 1 << i
    return n

def expand_gate(num_qubits, gate, to_qubits):
    N = 2**num_qubits
    num_to_qubits = len(to_qubits)
    num_not_to_qubits = num_qubits - num_to_qubits
    to_qubits_set = set(to_qubits)
    not_to_qubits_set = set(range(num_qubits)) - to_qubits_set
    not_to_qubits = sorted(not_to_qubits_set)
    #M = np.zeros((N, N), complex)
    M_data = []
    M_i_indexes = []
    M_j_indexes = []

    for base_state in range(2**num_not_to_qubits):
        base_bits = int_to_bits(base_state, num_not_to_qubits)
        i = np.empty(num_qubits, np.bool)
        for q, q_index in enumerate(not_to_qubits):
            i[q_index] = base_bits[q]
        j = i.copy()

        for i_star, j_star in product(range(2**num_to_qubits), range(2**num_to_qubits)):
            i_star_bits, j_star_bits = int_to_bits(i_star, 2**num_to_qubits), int_to_bits(j_star, 2**num_to_qubits)
            for M_index, G_index in zip((i,j), (i_star_bits, j_star_bits)):
                for G_bit_index, M_bit_index in enumerate(to_qubits):
                    M_index[M_bit_index] = G_index[G_bit_index]

            #M[bit_vector_to_int(i), bit_vector_to_int(j)] = gate[i_star, j_star]
            M_data.append(gate[i_star, j_star])
            M_i_indexes.append(bit_vector_to_int(i))
            M_j_indexes.append(bit_vector_to_int(j))

    return sparse.csr_matrix((M_data, (M_i_indexes, M_j_indexes)), shape=(N, N))

def apply_gate(register, gate):
    return register * gate

class QuantumComputer:
    def __init__(self, num_qubits, initial_state=0):
        self.num_qubits = num_qubits
        self.num_states = 2**num_qubits
        self.register = np.zeros(self.num_states, complex)
        self.used_gates = {}

        if type(initial_state) is int:
            self.register[initial_state] = 1
        elif type(initial_state) is list:
            assert len(initial_state) == num_qubits, "Initial state must be same length as the number of qubits"
            self.register[bit_vector_to_int(initial_state)] = 1

        self.gates = {}
        self.gates["Pauli_X"] = np.array([[0, 1],
                                          [1, 0]])
        self.gates["Puali_Y"] = np.array([[0, -1j],
                                          [1j, 0]])
        self.gates["Pauli_Z"] = np.array([[1, 0],
                                          [0, -1]])
        self.gates["Hadamard"] = 1/np.sqrt(2) * np.array([[1, 1],
                                                          [1, -1]])

    def define_gate(self, name, gate):
        self.gates[name] = gate

    def apply_gate(self, gate_name, to_qubits, controls=(), notted_controls=()):
        if type(to_qubits) is int:
            to_qubits = (to_qubits,)
        if (gate_name, to_qubits, controls) not in self.used_gates:
            gate = self.gates[gate_name]
            
            if controls:
                num_gate_qubits = len(to_qubits) + len(controls) # + len(notted_controls)
                controlled_gate = sparse.eye(2**num_gate_qubits, dtype=complex, format="lil")
                starting_index = 2**num_gate_qubits - 2**len(to_qubits)
                controlled_gate[starting_index:, starting_index:] = gate
                expanded_gate = expand_gate(self.num_qubits, controlled_gate, to_qubits + controls)
            else:
                expanded_gate = expand_gate(self.num_qubits, gate, to_qubits)

            self.used_gates[gate_name, to_qubits, controls] = expanded_gate

        self.register = self.used_gates[gate_name, to_qubits, controls] * self.register

    def measure(self):
        measured = np.square(np.abs(self.register))
        likely_state = np.argmax(measured)
        print("The most likely state is " + str(likely_state) + " (" + int_to_bit_string(likely_state, self.num_qubits) + ") with a probability of " + str(measured[likely_state]) + ".")
        return measured

qc = QuantumComputer(3)
qc.apply_gate("Pauli_X", 0)
#qc.apply_gate("Pauli_X", 1)
qc.apply_gate("Pauli_X", 1, controls=(2,))
#qc.apply_gate("Pauli_X", 2, controls=(3,))
qc.measure()

