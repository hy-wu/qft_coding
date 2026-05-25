# -*- coding: utf-8 -*-
"""
An educational, from-scratch Quantum Simulator using only NumPy.
Designed for physicists to see exactly what happens mathematically under the hood
of high-level frameworks like Qiskit.
"""

import numpy as np

# Set print options to format complex numbers beautifully
np.set_printoptions(precision=4, suppress=True)

class PhysicistQuantumSimulator:
    """
    A quantum simulator that represents the quantum state as a raw statevector
    in Hilbert space, and applies gates via direct matrix multiplication.
    """
    def __init__(self, num_qubits: int):
        self.num_qubits = num_qubits
        self.dim = 2 ** num_qubits
        
        # Initialize the state to |00...0>
        # In Hilbert space, this is a vector of size 2^N with a 1 at index 0 and 0 elsewhere.
        self.statevector = np.zeros(self.dim, dtype=complex)
        self.statevector[0] = 1.0
        
        # Standard Single-Qubit Gates in their 2x2 matrix representations
        self.I = np.array([[1.0, 0.0], 
                           [0.0, 1.0]], dtype=complex)
        
        self.H = np.array([[1.0, 1.0], 
                           [1.0, -1.0]], dtype=complex) / np.sqrt(2.0)
        
        self.X = np.array([[0.0, 1.0], 
                           [1.0, 0.0]], dtype=complex)
        
        self.Y = np.array([[0.0, -1.0j], 
                           [1.0j, 0.0]], dtype=complex)
        
        self.Z = np.array([[1.0, 0.0], 
                           [0.0, -1.0]], dtype=complex)

    def print_state(self, message: str = ""):
        """
        Prints the current statevector as a linear combination of Dirac ket states.
        """
        if message:
            print(f"\n--- {message} ---")
        
        terms = []
        for i in range(self.dim):
            amplitude = self.statevector[i]
            if np.abs(amplitude) > 1e-10:
                # Format binary representation of index i to match |q_{N-1}...q_0>
                binary = format(i, f'0{self.num_qubits}b')
                # Format the complex amplitude beautifully
                real = amplitude.real
                imag = amplitude.imag
                if np.abs(imag) < 1e-10:
                    amp_str = f"{real:.4f}"
                elif np.abs(real) < 1e-10:
                    amp_str = f"{imag:.4f}j"
                else:
                    amp_str = f"({real:.4f} + {imag:.4f}j)"
                
                terms.append(f"{amp_str}|{binary}>")
        
        state_str = " + ".join(terms)
        print(f"Statevector: {self.statevector}")
        print(f"Ket Notation:  |ψ> = {state_str}")

    def apply_single_qubit_gate(self, gate_matrix: np.ndarray, target: int):
        """
        Applies a 2x2 single-qubit gate to a specific target qubit.
        Under the hood, this constructs a 2^N x 2^N matrix using Kronecker tensor products.
        Qiskit convention (LSB-first): the state is ordered as |q_{N-1} ... q_1 q_0>
        So Kronecker product runs from N-1 down to 0:
        G = I^(N-1) ⊗ ... ⊗ Gate^(target) ⊗ ... ⊗ I^(0)
        """
        operator = np.array([[1.0]], dtype=complex)
        for j in reversed(range(self.num_qubits)):
            if j == target:
                operator = np.kron(operator, gate_matrix)
            else:
                operator = np.kron(operator, self.I)
        
        # Apply the full operator via matrix multiplication: |ψ_new> = G |ψ_old>
        self.statevector = operator @ self.statevector

    def apply_cnot(self, control: int, target: int):
        """
        Applies a Controlled-NOT (CNOT) gate between control and target qubits.
        In Hilbert space, CNOT is a 2^N x 2^N permutation matrix where
        the state of the target qubit is flipped if and only if the control qubit is 1.
        """
        cnot_matrix = np.zeros((self.dim, self.dim), dtype=complex)
        
        for x in range(self.dim):
            # Read the control bit of the basis state x
            control_bit = (x >> control) & 1
            if control_bit == 1:
                # Flip the target bit of x to get the new basis state y
                y = x ^ (1 << target)
            else:
                y = x
            cnot_matrix[y, x] = 1.0
            
        # Apply CNOT: |ψ_new> = CNOT |ψ_old>
        self.statevector = cnot_matrix @ self.statevector

    def measure(self, shots: int = 1024):
        """
        Simulates measurement of all qubits by sampling from the probability distribution
        dictated by the Born Rule: P(x) = |<x|ψ>|^2.
        """
        # Calculate probabilities from amplitudes (Born Rule)
        probabilities = np.abs(self.statevector) ** 2
        
        # Randomly sample based on these probabilities
        outcomes = np.random.choice(self.dim, size=shots, p=probabilities)
        
        # Count the frequency of each outcome
        counts = {}
        for outcome in outcomes:
            binary = format(outcome, f'0{self.num_qubits}b')
            counts[binary] = counts.get(binary, 0) + 1
            
        return counts


if __name__ == "__main__":
    N = 12
    print("======================================================================")
    print(f"      PHYSISIST'S FROM-SCRATCH QUANTUM SIMULATOR ({N} QUBITS)           ")
    print("======================================================================")
    
    # 1. Initialize our N-qubit system
    sim = PhysicistQuantumSimulator(num_qubits=N)
    sim.print_state(f"Initial State (Ground State |{'0' * N}>)")
    
    # 2. Apply Hadamard to qubit 0 (qc.h(0))
    # H gate maps |0> -> (|0> + |1>)/sqrt(2).
    # Since qubit 1 remains in |0>, the state becomes:
    # (|0>_1 ⊗ |0>_0 + |0>_1 ⊗ |1>_0)/sqrt(2) = (|00> + |01>)/sqrt(2)
    sim.apply_single_qubit_gate(sim.H, target=0)
    sim.print_state(f"After applying H on qubit 0 (Superposition) |{'0' * (N-1)}1>")
    
    # 3. Apply CNOT with control = 0, target = 1 (qc.cx(0, 1))
    # CNOT is the entangling operation:
    # If qubit 0 (control) is 0 -> qubit 1 (target) remains 0. (|00> remains |00>)
    # If qubit 0 (control) is 1 -> qubit 1 (target) flips from 0 to 1. (|01> becomes |11>)
    # This creates the Bell State: (|00> + |11>)/sqrt(2)
    sim.apply_cnot(control=0, target=1)
    sim.print_state("After CNOT(control=0, target=1) (Entanglement/Bell State)")
    
    # Let's verify mathematically why this state is entangled:
    # A state is entangled if it CANNOT be written as a tensor product of two separate qubits:
    # |ψ> != (a|0> + b|1>) ⊗ (c|0> + d|1>) = ac|00> + ad|01> + bc|10> + bd|11>
    # For Bell state: ac = 1/sqrt(2), ad = 0, bc = 0, bd = 1/sqrt(2).
    # ad = 0 implies either a=0 or d=0.
    # If a=0, then ac=0, which contradicts ac = 1/sqrt(2).
    # Thus, this state is mathematically inseparable — this is a purely quantum phenomenon!
    
    # 4. Simulate Measurement
    print("\n--- Simulating 1024 measurements (Born Rule: P(x) = |<x|ψ>|^2) ---")
    counts = sim.measure(shots=1024)
    print(f"Measurement Counts: {counts}")
