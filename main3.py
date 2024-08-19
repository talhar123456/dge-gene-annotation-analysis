import numpy as np
import logomaker
import matplotlib.pyplot as plt

class PSSM:
    def __init__(self, sequences):
        self.sequences = sequences
        self.alphabet = ['A', 'C', 'G', 'T']
        self.l = len(sequences)
        self.seq_length = len(sequences[0])

    def get_FrequencyMatrix(self):
        freq_matrix = np.zeros((self.seq_length, len(self.alphabet)))
        for seq in self.sequences:
            for j, nucleotide in enumerate(seq):
                i = self.alphabet.index(nucleotide)
                freq_matrix[j, i] += 1
        freq_matrix /= self.l
        return freq_matrix

    def get_CorrectedFrequencyMatrix(self, pseudo_weights):
        freq_matrix = self.get_FrequencyMatrix()
        corrected_freq_matrix = np.zeros_like(freq_matrix)
        k = sum(pseudo_weights.values())
        for i, nucleotide in enumerate(self.alphabet):
            corrected_freq_matrix[:, i] = (freq_matrix[:, i] + pseudo_weights[nucleotide]) / (self.l + k)
        return corrected_freq_matrix

    def get_ScoringMatrix(self, corrected_freq_matrix, pseudo_weights):
        scoring_matrix = np.zeros_like(corrected_freq_matrix)
        for i, nucleotide in enumerate(self.alphabet):
            scoring_matrix[:, i] = corrected_freq_matrix[:, i] / pseudo_weights[nucleotide]
        return scoring_matrix

    def plot_SequenceLogo(self, matrix, title):
        matrix_df = pd.DataFrame(matrix, columns=self.alphabet)
        logomaker.Logo(matrix_df, shade_below=.5, fade_below=.5, font_name='Arial Rounded MT Bold')
        plt.title(title)
        plt.xlabel('Position')
        plt.ylabel('Frequency')
        plt.show()

import pandas as pd

sequences = [
    "TCACACGTGGGA",
    "GGCCACGTGCAG",
    "TGACACGTGGGA",
    "CAGCACGTGGGG",
    "TACCACGTGCGA",
    "ACGCACGTTGGT",
    "CAGCACGTTTTC",
    "TAGCACGTTTTC"
]

pseudo_weights = {'A': 0.3, 'T': 0.3, 'G': 0.2, 'C': 0.2}

pssm = PSSM(sequences)
frequency_matrix = pssm.get_FrequencyMatrix()
print("Frequency Matrix:\n", frequency_matrix)

corrected_freq_matrix = pssm.get_CorrectedFrequencyMatrix(pseudo_weights)
print("Corrected Frequency Matrix:\n", corrected_freq_matrix)

scoring_matrix = pssm.get_ScoringMatrix(corrected_freq_matrix, pseudo_weights)
print("Scoring Matrix:\n", scoring_matrix)

pssm.plot_SequenceLogo(frequency_matrix, 'Frequency Matrix')
pssm.plot_SequenceLogo(corrected_freq_matrix, 'Corrected Frequency Matrix')
pssm.plot_SequenceLogo(scoring_matrix, 'Scoring Matrix')
