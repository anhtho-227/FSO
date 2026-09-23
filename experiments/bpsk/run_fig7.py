"""
run_fig7.py -- Reconstruct Figure 7: Analytical vs Simulated (bit stream simulation) (Section IX.E)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import TABLE_I, Pe1_closed_form_eq39
from data_generate.scripts.generate_channel_samples import generate_bpsk_bitstream

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")

def main():
    sigma = 0.5
    gammas = np.linspace(1, 6, 8)
    
    fig, ax = plt.subplots(figsize=(7, 5))
    for idx, (regime, (a, b)) in enumerate(TABLE_I.items()):
        pe_ana = [Pe1_closed_form_eq39(g, a, b, sigma) for g in gammas]
        pe_sim = []
        for g in gammas:
            bits, I1, eps1 = generate_bpsk_bitstream(a, b, sigma, 300_000, seed=42+idx)
            psi = np.where(bits == 0, g, -g)
            pe_cond = np.clip(0.25 * np.exp(-4 * (np.sqrt(I1) * psi + eps1) ** 2), 0, 1)
            rng = np.random.default_rng(99)
            errors = rng.random(len(bits)) < pe_cond
            pe_sim.append(errors.mean())
            
        ax.semilogy(gammas, pe_ana, '-', label=f'{regime} (analytical)')
        ax.semilogy(gammas, pe_sim, 's', markerfacecolor='none', label=f'{regime} (simulated)')
        
    ax.set_xlabel('|γ|')
    ax.set_ylabel('Helstrom Error Probability (log)')
    ax.set_title('Figure 7 — Analytical vs Simulated (3×10^5 bits/point)')
    ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(fontsize=7)
    
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "Fig7_reconstructed.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Saved: {path}")

if __name__ == "__main__":
    main()
