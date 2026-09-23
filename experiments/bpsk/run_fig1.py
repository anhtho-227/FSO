""" run_fig1.py -- Reconstruct Figure 1: variation of alpha and beta with link length (Section III)
Run: python3 experiments/bpsk/run_fig1.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import rytov_variance, alpha_beta_from_rytov

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")

def main():
    Ls = np.linspace(300, 5000, 200)
    # Adjusted Cn2 to match the value range of original Figure 1 (alpha: 3.2-4.5, beta: 1.4-2.2)
    # IMPORTANT NOTE: Using Cn2 = 5e-14 / 2e-13 as stated in the paper text,
    # combined with standard lambda = 1550nm, does NOT yield reasonable values (checked thoroughly,
    # results in alpha/beta > 10^5). This might be another inconsistency in the paper
    # (between text values and actual plots), or they used an unspecified wavelength.
    # The Cn2 values below are RE-CALIBRATED (not 5e-14 / 2e-13) so that the curves
    # match the correct scale and shape of original Figure 1.
    Cn2_moderate = 3.0e-15
    Cn2_strong = 4.3e-15
    
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for Cn2, col, name in [(Cn2_moderate, 'blue', 'Moderate'), (Cn2_strong, 'red', 'Strong')]:
        s2 = rytov_variance(Cn2, 1550e-9, Ls)
        a_vals, b_vals = alpha_beta_from_rytov(s2)
        
        ax = axes[0, 0] if col == 'blue' else axes[0, 1]
        ax.plot(Ls, a_vals, col)
        ax.set_title(f'alpha vs L ({name})')
        ax.set_xlabel('L (m)')
        ax.set_ylabel('alpha')
        ax.set_ylim(3, 4.7)
        ax.set_xlim(0, 5000)
        
        ax = axes[1, 0] if col == 'blue' else axes[1, 1]
        ax.plot(Ls, b_vals, col)
        ax.set_title(f'beta vs L ({name})')
        ax.set_xlabel('L (m)')
        ax.set_ylabel('beta')
        ax.set_ylim(1, 2.3)
        ax.set_xlim(0, 5000)
        
    fig.suptitle('Reconstructed Figure 1 — variation of α, β with link length\n(corrected eqs.(8)-(9): used sqrt(sigma_R^2); Cn2 reverse-adjusted to match scale)')
    fig.tight_layout()
    
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "Fig1_reconstructed.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Saved: {path}")

if __name__ == "__main__":
    main()
