import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import random

# Mock data generation for visualization
def generate_mock_data(n=200):
    data = []
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    
    for _ in range(n):
        # Generate random sequence
        length = random.randint(10, 20)
        seq = "".join(random.choice(amino_acids) for _ in range(length))
        
        # Calculate Hydrophobicity (GRAVY)
        analysis = ProteinAnalysis(seq)
        gravy = analysis.gravy()
        
        # Mock Affinity (correlated slightly with hydrophobicity for the "trap" effect, or just random)
        # Often high affinity binders can be hydrophobic.
        affinity = random.normalvariate(0.5, 0.2)
        affinity = max(0, min(1, affinity)) # Clamp to 0-1
        
        data.append({'sequence': seq, 'gravy': gravy, 'affinity': affinity})
    
    return pd.DataFrame(data)

def plot_selection():
    df = generate_mock_data(300)
    
    # Thresholds
    gravy_threshold = 0.5 # Max allowed hydrophobicity
    
    # Classify
    # "Green Box": Low Hydrophobicity (<= 0.5) AND High Affinity (e.g., > 0.6)
    # But the filter only strictly cuts on GRAVY.
    # The "Pitch" is: We want High Affinity but Low Hydrophobicity.
    
    df['status'] = 'Rejected'
    df.loc[df['gravy'] <= gravy_threshold, 'status'] = 'Passed'
    
    # Create plot
    plt.figure(figsize=(10, 6))
    
    # Scatter plot
    passed = df[df['status'] == 'Passed']
    rejected = df[df['status'] == 'Rejected']
    
    plt.scatter(rejected['gravy'], rejected['affinity'], c='red', alpha=0.5, label='Rejected (High Exhaustion Risk)')
    plt.scatter(passed['gravy'], passed['affinity'], c='green', alpha=0.6, label='Passed (TME Stable)')
    
    # Highlight the "Green Box" (Target Region)
    # Low Hydrophobicity (< 0.5), High Affinity (> 0.5)
    plt.axvline(x=gravy_threshold, color='black', linestyle='--', label='Max Hydrophobicity Threshold')
    
    # Add a shaded green region for the "Ideal" candidates
    # Assuming ideal is GRAVY < 0.5 and Affinity > 0.5
    plt.fill_betweenx([0.5, 1.0], -2.0, gravy_threshold, color='green', alpha=0.1, label='Target Region')

    plt.title('TCR Candidate Selection: TME Survival vs. Affinity')
    plt.xlabel('Tonic Signaling Risk (Hydrophobicity / GRAVY)')
    plt.ylabel('Predicted Affinity')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Save plot
    plt.savefig('results/selection_plot.png')
    print("Plot saved to results/selection_plot.png")

if __name__ == "__main__":
    plot_selection()
