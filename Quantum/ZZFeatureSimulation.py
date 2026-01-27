import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from qiskit.circuit.library import ZZFeatureMap
from qiskit_machine_learning.kernels import FidelityQuantumKernel
from sklearn.preprocessing import MinMaxScaler
from tqdm import tqdm


class ZZFeatureSimulation:


     def __init__(self,
                  data: pd.DataFrame,
                  map_range: float)->None:
         self.data = data
         self.num_qubits = len(data.columns)
         self.encoding_range = map_range

         scaler = MinMaxScaler(feature_range=(0, 1))
         X_scaled_01 = scaler.fit_transform(self.data)

         # Ora espandiamo al range quantistico [0, we]
         X_quantum = X_scaled_01 * self.encoding_range

         feature_map = ZZFeatureMap(feature_dimension=self.num_qubits, reps=1, entanglement='linear')

         kernel = FidelityQuantumKernel(feature_map=feature_map)

         # print("Calcolo della matrice Kernel in corso... ")
         # kernel_matrix = kernel.evaluate(x_vec=X.values)
         n_pazienti = len(self.data)
         self.kernel_matrix = np.zeros((n_pazienti, n_pazienti))
         data_values = X_quantum

         print(f"Inizio calcolo Kernel Quantistico per {n_pazienti} pazienti.")
         print("Calcolando solo la metà superiore della matrice per simmetria...")

         for i in tqdm(range(n_pazienti), desc="Avanzamento"):
             for j in range(i, n_pazienti):  # Parte da i: calcoliamo solo il triangolo superiore

                 if i == j:
                     # La diagonale è sempre 1.0 per definizione
                     self.kernel_matrix[i, j] = 1.0
                 else:
                     # Calcoliamo la somiglianza tra Paziente i e Paziente j
                     # Passiamo solo due righe al kernel
                     paziente_A = data_values[i].reshape(1, -1)
                     paziente_B = data_values[j].reshape(1, -1)

                     # evaluate restituisce una matrice 1x1, prendiamo il valore scalare
                     valore = kernel.evaluate(x_vec=paziente_A, y_vec=paziente_B)[0, 0]

                     # Riempiamo la matrice (Simmetria: A vs B è uguale a B vs A)
                     self.kernel_matrix[i, j] = valore
                     self.kernel_matrix[j, i] = valore

         print("\nCalcolo completato!")
         print(f"Dimensioni Matrice: {self.kernel_matrix.shape}")

     def get_kernel_matrix(self):
         return self.kernel_matrix

     def get_encoding_range(self):
         return self.encoding_range

     def spectral_clustering(self, n_clusters):
         from sklearn.cluster import SpectralClustering
         self.clustering = SpectralClustering(n_clusters=n_clusters, affinity='precomputed')
         labels = self.clustering.fit_predict(self.kernel_matrix)
         return labels

     def plot_kernel_matrix(self):
         import seaborn as sns
         sorted_indices = np.argsort(self.spectral_clustering(3))

         # 2. Riordina la Matrice Kernel (sia righe che colonne)
         sorted_kernel = self.kernel_matrix[sorted_indices][:, sorted_indices]

         viz_kernel = sorted_kernel.copy()

         # 2. "Cancelliamo" la diagonale
         # Questo costringe la scala colori a concentrarsi sulle differenze reali tra pazienti
         np.fill_diagonal(viz_kernel, np.nan)

         plt.figure(figsize=(10, 8))
         # 'robust=True' aiuta a ignorare eventuali outlier che rovinano la scala
         sns.heatmap(viz_kernel, cmap='viridis', xticklabels=False, yticklabels=False, robust=True)
         plt.title("Matrice Kernel (Diagonale Rimossa per Contrasto)")
         plt.show()