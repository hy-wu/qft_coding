from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFT
import numpy as np

def build_qft_circuit(num_qubits: int) -> QuantumCircuit:
    qc = QuantumCircuit(num_qubits)
    for i in range(num_qubits):
        qc.h(i)
        for j in range(i + 1, num_qubits):
            k = j - i + 1
            theta = 2 * np.pi / (2**k)
            qc.cp(theta, j, i)
    for i in range(num_qubits // 2):
        qc.swap(i, num_qubits - 1 - i)
    return qc

# 1. Custom QFT
qc_custom = build_qft_circuit(4)
print("=== CUSTOM QFT CIRCUIT ===")
print(qc_custom.draw(output='text'))

# 2. Qiskit Library QFT
qc_lib = QFT(num_qubits=4)
print("\n=== QISKIT LIBRARY QFT CIRCUIT ===")
print(qc_lib.decompose().draw(output='text'))
