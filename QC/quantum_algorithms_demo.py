# -*- coding: utf-8 -*-
"""
Demonstration of core Quantum Algorithms:
1. Quantum Fourier Transform (QFT) on 4 Qubits (Frequency Analysis).
2. Grover's Search Algorithm on 3 Qubits (Amplitude Amplification).
Uses Qiskit and AerSimulator to demonstrate physical step-by-step outputs.
"""

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator

def get_statevector(qc: QuantumCircuit) -> np.ndarray:
    """Helper to extract the exact statevector from a circuit."""
    qc_sv = qc.copy()
    qc_sv.remove_final_measurements() # Remove measurements to get pure state
    qc_sv.save_statevector()
    
    sim = AerSimulator(method='statevector')
    t_qc = transpile(qc_sv, sim)
    result = sim.run(t_qc).result()
    return np.array(result.get_statevector())

# ======================================================================
# 1. QUANTUM FOURIER TRANSFORM (QFT) DEMO
# ======================================================================
def build_qft_circuit(num_qubits: int) -> QuantumCircuit:
    """Builds the standard QFT circuit for a given number of qubits."""
    qc = QuantumCircuit(num_qubits)
    
    for i in reversed(range(num_qubits)):
        # 1. Apply Hadamard on the active qubit
        qc.h(i)
        # 2. Apply controlled phase rotations
        for j in reversed(range(i)):
            k = i - j + 1
            theta = 2 * np.pi / (2**k)
            qc.cp(theta, j, i)
            
    # 3. Swap qubits at the end to match standard binary ordering
    for i in range(num_qubits // 2):
        qc.swap(i, num_qubits - 1 - i)
        
    return qc

def run_qft_demo():
    print("======================================================================")
    print("      1. QUANTUM FOURIER TRANSFORM (QFT) ON 4 QUBITS                  ")
    print("======================================================================")
    
    N = 4
    # Create a periodic state in the computational basis.
    # Let's prepare a superposition that repeats with a certain period.
    # Specifically, a state of period 4: |ψ> = (|0000> + |0100> + |1000> + |1100>) / 2
    # In decimal indices, this is a superposition of states: 0, 4, 8, 12.
    qc = QuantumCircuit(N)
    
    # Prepare the periodic state:
    # Qubits 0 and 1 are in |0>
    # Qubits 2 and 3 are put in superposition using Hadamard:
    # This creates a periodic signal across the 16 computational states!
    qc.h(2)
    qc.h(3)
    
    initial_state = get_statevector(qc)
    print("Initial periodic state |ψ_in> (Indices with non-zero amplitudes):")
    for idx, amp in enumerate(initial_state):
        if np.abs(amp) > 1e-5:
            print(f"  |{idx:04b}> (dec {idx:2d}): Amplitude = {amp.real:.4f} + {amp.imag:.4f}j")
            
    # Append the QFT circuit
    qft_qc = build_qft_circuit(N)
    qc.compose(qft_qc, inplace=True)
    
    final_state = get_statevector(qc)
    print("\nAfter QFT |ψ_out> (Fourier transformed frequencies):")
    # A period 4 signal in a 16-dimensional space should yield delta peaks at 
    # frequency indices that are multiples of 16 / 4 = 4.
    # i.e., at dec indices 0, 4, 8, 12.
    for idx, amp in enumerate(final_state):
        if np.abs(amp) > 1e-5:
            print(f"  |{idx:04b}> (dec {idx:2d}): Amplitude = {amp.real:.4f} + {amp.imag:.4f}j")
            
    print("\n[PHYSICAL ANALYSIS] Notice that the periodic state input in the spatial domain")
    print("has been mapped directly to clean delta-like frequency peaks in the Fourier domain.")
    print("This is the exact quantum equivalent of the classical Discrete Fourier Transform (DFT)!")


# ======================================================================
# 2. GROVER'S SEARCH ALGORITHM DEMO
# ======================================================================
def add_grover_oracle(qc: QuantumCircuit, target_state: str):
    """
    Applies the Phase Oracle for a target state on 3 qubits.
    Flips the sign of the target state: |x> -> -|x> if x == target, else |x>.
    For 3 qubits, we use a Controlled-Controlled-Z (CCZ) gate.
    """
    num_qubits = qc.num_qubits
    # Apply X gates to match target state bitstring '1's and '0's
    # Qiskit is LSB-first, so target_state is read right-to-left
    for i, bit in enumerate(reversed(target_state)):
        if bit == '0':
            qc.x(i)
            
    # Apply Controlled-Controlled-Z (CCZ) on 3 qubits
    # CCZ = H on target -> CCNOT (Toffoli) -> H on target
    qc.h(2)
    qc.ccx(0, 1, 2)
    qc.h(2)
    
    # Undo the X gates to restore computational basis
    for i, bit in enumerate(reversed(target_state)):
        if bit == '0':
            qc.x(i)

def add_grover_diffuser(qc: QuantumCircuit):
    """
    Applies the Grover Diffuser (inversion about the mean) on 3 qubits.
    U_s = 2|s><s| - I = H^⊗3 (2|0><0| - I) H^⊗3
    """
    num_qubits = qc.num_qubits
    for i in range(num_qubits):
        qc.h(i)
        qc.x(i)
        
    # Apply CCZ (Phase flip on |000> state)
    qc.h(2)
    qc.ccx(0, 1, 2)
    qc.h(2)
    
    for i in range(num_qubits):
        qc.x(i)
        qc.h(i)

def run_grover_demo():
    print("\n======================================================================")
    print("      2. GROVER'S SEARCH ALGORITHM ON 3 QUBITS (TARGET: |101>)        ")
    print("======================================================================")
    
    N = 3
    target_state = "101" # Qubit 2 = 1, Qubit 1 = 0, Qubit 0 = 1 (decimal index 5)
    
    # Step 1: Initialize in a uniform superposition of all 2^3 = 8 states
    qc = QuantumCircuit(N)
    for i in range(N):
        qc.h(i)
        
    state = get_statevector(qc)
    print(f"Step 0: Initial Uniform Superposition (Probability of finding each state = 12.5%):")
    print(f"  P({target_state}) = {np.abs(state[5])**2 * 100:.2f}%")
    
    # Optimal number of steps for 3 qubits (N=8 states): T ≈ (pi/4)*sqrt(8) ≈ 2.22 -> 2 steps
    for step in range(1, 3):
        # Apply Oracle (flips phase of target state)
        add_grover_oracle(qc, target_state)
        # Apply Diffuser (inverts amplitudes about the mean)
        add_grover_diffuser(qc)
        
        state = get_statevector(qc)
        print(f"\nStep {step}: After {step} Grover Iteration(s):")
        print(f"  Amplitudes of all 8 states:")
        for idx, amp in enumerate(state):
            binary = f"{idx:03b}"
            is_target = " <-- TARGET" if binary == target_state else ""
            prob = np.abs(amp)**2 * 100
            print(f"    |{binary}> (dec {idx}): Amplitude = {amp.real:.4f}, Prob = {prob:.2f}% {is_target}")
            
    print("\n[PHYSICAL ANALYSIS] Notice how the probability of our target state |101>")
    print("amplified from 12.5% (uniform) to 94.53% in just 2 steps!")
    print("Grover's algorithm achieves a quadratic quantum speedup by turning unstructured")
    print("search into a coherent wave interference process in Hilbert space.")

if __name__ == "__main__":
    run_qft_demo()
    run_grover_demo()
