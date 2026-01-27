import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import joblib
from qiskit import generate_preset_pass_manager, transpile

from sklearn.preprocessing import MinMaxScaler
from sklearn.cluster import SpectralClustering
from sklearn.metrics.pairwise import rbf_kernel
import os
from dotenv import load_dotenv

# Qiskit Imports
from qiskit_ibm_runtime import QiskitRuntimeService, Session
from qiskit_ibm_runtime import SamplerV2 as RuntimeSampler
from qiskit_ibm_runtime import EstimatorV2 as RuntimeEstimator
from qiskit.circuit.library import ZZFeatureMap, zz_feature_map
from qiskit_algorithms.state_fidelities import ComputeUncompute
from qiskit_machine_learning.kernels import FidelityQuantumKernel
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator

# -----------------------------------------------------------
# 1. CONFIGURAZIONE CLOUD E BACKEND
# -----------------------------------------------------------


load_dotenv()


api_key = os.getenv("IBM_CLOUD_API_KEY")
crn = os.getenv("IBM_CLOUD_CRN")


service = QiskitRuntimeService(
    channel="ibm_quantum_platform",
    token=api_key,
    instance=crn
)


backend = service.least_busy(operational=True, simulator=False, min_num_qubits=7)

print(f"🚀 Connesso al backend: {backend.name}")


df = pd.read_csv('quantum-alz.csv', sep=";", decimal=',')

X = df[['Cope_Avoidant', 'Cope_Active', 'BDI', 'SES', 'EHI', 'NEO_CON', 'NEO_AGR', 'AUDIT', 'Cope_Emotional', 'CVLT_9', 'N_CVLT_1', 'N_CVLT_2', 'Dim_1', 'Dim_4']]

num_qubits = len(X.columns)
scaler = MinMaxScaler(feature_range=(0, 1))
X_normalized = scaler.fit_transform(X)

# ==========================================
# 0. PASSAGGIO CRUCIALE: TRANSPILAZIONE
# ==========================================
print("🔧 Transpilazione circuiti per l'hardware specifico...")

# Creiamo la Feature Map astratta
abstract_circuit = zz_feature_map(feature_dimension=num_qubits, reps=2, entanglement='linear')


pass_manager = generate_preset_pass_manager(
    optimization_level=3, backend=backend
)

# Transpiliamo il circuito ORA (diventa un circuito ISA)
isa_feature_map = pass_manager.run(abstract_circuit)

print("✅ Circuito convertito in ISA (Instruction Set Architecture)!")

# ==========================================
# 1. FQK (Fidelity Quantum Kernel)
# ==========================================
'''print("\n--- INIZIO CALCOLO FQK ---")

SCALE_FQK = 0.5
X_fqk = X_normalized * SCALE_FQK

# Configura Sampler in Job Mode
sampler = RuntimeSampler(mode=backend)
sampler.options.default_shots = 4096

# Collega il Sampler all'algoritmo
fidelity = ComputeUncompute(sampler=sampler)

# IMPORTANTE: Passiamo 'isa_feature_map', non quello astratto!
fqk_kernel = FidelityQuantumKernel(feature_map=isa_feature_map, fidelity=fidelity)

print(f"⏳ Invio job FQK a {backend.name}...")
# Nota: FidelityQuantumKernel gestisce internamente i parametri,
# ma il circuito base deve essere già transpilato.
matrix_fqk = fqk_kernel.evaluate(x_vec=X_fqk)

print("✅ Matrice FQK calcolata!")

# D. Clustering FQK
print("   Eseguo Spectral Clustering su FQK...")
sc_fqk = SpectralClustering(n_clusters=3, affinity='precomputed', assign_labels='discretize', random_state=42)
labels_fqk = sc_fqk.fit_predict(matrix_fqk)
'''
# ==========================================
# 2. PQK (Projected Quantum Kernel)
# ==========================================
print("\n--- INIZIO CALCOLO PQK ---")

SCALE_PQK = 2 * np.pi
X_pqk = X_normalized * SCALE_PQK

# Definiamo gli osservabili astratti
observables = [SparsePauliOp.from_list([("I" * i + "Z" + "I" * (num_qubits - 1 - i), 1)])
               for i in range(num_qubits)]

# CRUCIALE PER PQK: Anche gli osservabili devono essere adattati al layout del circuito transpilato
print("🔧 Adattamento osservabili al layout fisico...")
isa_observables = [obs.apply_layout(isa_feature_map.layout) for obs in observables]


def run_pqk_job_mode(data, transpiled_map, isa_obs, backend_obj):
    # 1. Prepariamo i PUBs usando il circuito GIÀ transpilato
    pubs = []
    for x in data:
        # Assegniamo i parametri al circuito transpilato
        bound_circ = transpiled_map.assign_parameters(x)
        pubs.append((bound_circ, isa_obs))

    print(f"⏳ Invio job PQK ({len(data)} circuiti) a {backend_obj.name}...")

    # 2. Estimator V2
    estimator = RuntimeEstimator(mode=backend_obj)
    estimator.options.default_shots = 4096

    # 3. Run
    job = estimator.run(pubs)
    print(f"   Job ID: {job.job_id()}")

    results = job.result()
    features = [res.data.evs for res in results]
    return np.array(features)


# Esecuzione usando feature map e osservabili ISA
features_pqk = run_pqk_job_mode(X_pqk, isa_feature_map, isa_observables, backend)

# Calcolo Kernel Classico
gamma_val = 1.0 / num_qubits
matrix_pqk = rbf_kernel(features_pqk, gamma=gamma_val)

print("✅ Matrice PQK calcolata!")

'''
# -----------------------------------------------------------
# 2. APPROCCIO A: FIDELITY QUANTUM KERNEL (FQK)
#    Parametro chiave: Scaling = 0.5
# -----------------------------------------------------------
print("\n--- INIZIO CALCOLO FQK (Scaling 0.5) ---")

# A. Scaling specifico per FQK
SCALE_FQK = 0.5
X_fqk = X_normalized * SCALE_FQK

# B. Feature Map
feature_map = ZZFeatureMap(feature_dimension=num_qubits, reps=2, entanglement='linear')

# C. Calcolo su Cloud con Sessione
# Usiamo SamplerV2 perché FQK misura la probabilità di misurare '000...0'

sampler = RuntimeSampler(mode=backend)
sampler.options.default_shots = 4096  # Precisione statistica

# Collega il Sampler del cloud all'algoritmo di Fedeltà
fidelity = ComputeUncompute(sampler=sampler)

# Crea il Kernel
fqk_kernel = FidelityQuantumKernel(feature_map=feature_map, fidelity=fidelity)

print("⏳ Invio job FQK al Cloud... (Questo calcola la matrice N x N)")
# evaluate calcola la matrice di similarità completa
matrix_fqk = fqk_kernel.evaluate(x_vec=X_fqk)

print("Matrice FQK calcolata!")

# D. Clustering FQK
print("   Eseguo Spectral Clustering su FQK...")
sc_fqk = SpectralClustering(n_clusters=3, affinity='precomputed', assign_labels='discretize', random_state=42)
labels_fqk = sc_fqk.fit_predict(matrix_fqk)

# -----------------------------------------------------------
# 3. APPROCCIO B: PROJECTED QUANTUM KERNEL (PQK)
#    Parametro chiave: Scaling = 2*pi
# -----------------------------------------------------------
print("\n--- INIZIO CALCOLO PQK (Scaling 2*pi) ---")

# A. Scaling specifico per PQK (Necessita rotazioni complete)
SCALE_PQK = 2 * np.pi
X_pqk = X_normalized * SCALE_PQK

# B. Osservabili (Magnetizzazione Z su tutti i qubit)
observables = [SparsePauliOp.from_list([("I" * i + "Z" + "I" * (num_qubits - 1 - i), 1)])
               for i in range(num_qubits)]


# C. Funzione Batching per EstimatorV2
def run_pqk_on_cloud(data, feature_map, obs, backend_service):
    # 1. Prepara i circuiti (PUBs)
    pubs = []
    for x in data:
        bound_circ = feature_map.assign_parameters(x)
        pubs.append((bound_circ, obs))

    # 2. Esegui in sessione
    
    estimator = RuntimeEstimator(mode=backend)
    estimator.options.default_shots = 4096

    print(f"⏳ Invio job PQK batch ({len(data)} circuiti) al Cloud...")
    job = estimator.run(pubs)
    results = job.result()

    # 3. Estrai dati (Expectation Values)
    features = [res.data.evs for res in results]
    return np.array(features)


# D. Esecuzione
features_pqk = run_pqk_on_cloud(X_pqk, feature_map, observables, backend)

# E. Calcolo Kernel Classico (RBF) sulle feature quantistiche
gamma_val = 1.0 / np.median(pd.DataFrame(features_pqk).var())  # Gamma euristica o 1/n_features
matrix_pqk = rbf_kernel(features_pqk, gamma=1.0 / num_qubits)

print("Matrice PQK calcolata!")'''

# F. Clustering PQK
print("   Eseguo Spectral Clustering su PQK...")
sc_pqk = SpectralClustering(n_clusters=3, affinity='precomputed', assign_labels='discretize', random_state=42)
labels_pqk = sc_pqk.fit_predict(matrix_pqk)

# -----------------------------------------------------------
# 4. VISUALIZZAZIONE COMPARATIVA
# -----------------------------------------------------------
plt.figure(figsize=(16, 7))

# Plot FQK
#sns.heatmap(matrix_fqk, cmap='viridis', ax=ax1)
#ax1.set_title(f"FQK (Scaling=0.5)\nBackend: {backend.name}")

# Plot PQK
sns.heatmap(matrix_pqk, cmap='viridis')
plt.title(f"PQK (Scaling=2$\\pi$)\nBackend: {backend.name}")

plt.show()


# 1. Salvataggio matrici
#joblib.dump(matrix_fqk, 'matrix_fqk.joblib')
joblib.dump(matrix_pqk, 'data/real_quantum/matrix_pqk.joblib')

# 2. Salvataggio cluster
#joblib.dump(sc_fqk, 'model_sc_fqk.joblib')
joblib.dump(sc_pqk, 'data/real_quantum/model_sc_pqk.joblib')

# 3. CRUCIALE: Salva i risultati già pronti (Labels e Indici)
# Creiamo un piccolo dizionario o dataframe per non perdere l'allineamento coi pazienti
results_data = {
    #'labels_fqk': labels_fqk,
    'labels_pqk': labels_pqk,
    'patient_indices': X.index  # Salviamo l'indice originale per ricollegarlo ai dati clinici
}
joblib.dump(results_data, 'data/real_quantum/cluster_results.joblib')

print("✅ Tutto salvato! Ora puoi chiudere questo script e spegnere il backend.")