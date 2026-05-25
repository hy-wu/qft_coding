import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit.library import QFT
from qiskit_aer import AerSimulator

def get_statevector(qc: QuantumCircuit) -> np.ndarray:
    qc_sv = qc.copy()
    qc_sv.remove_final_measurements()
    qc_sv.save_statevector()
    
    sim = AerSimulator(method='statevector')
    t_qc = transpile(qc_sv, sim)
    result = sim.run(t_qc).result()
    return np.array(result.get_statevector())

# Initialize periodic state: |++00>
# Superposition of 0, 4, 8, 12
qc_lib = QuantumCircuit(4)
qc_lib.h(2)
qc_lib.h(3)

# Apply Qiskit's built-in QFT
qc_lib.compose(QFT(num_qubits=4), inplace=True)

state_lib = get_statevector(qc_lib)
print("=== Qiskit Library QFT Output ===")
for idx, amp in enumerate(state_lib):
    if np.abs(amp) > 1e-5:
        print(f"  |{idx:04b}> (dec {idx:2d}): Amplitude = {amp.real:.4f} + {amp.imag:.4f}j")
