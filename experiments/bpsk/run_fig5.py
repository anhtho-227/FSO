"""
run_fig5.py -- Tái lập Hình 5: 3 chế độ turbulence x 2 mức nhiễu (Section IX.C)
Chạy: python3 experiments/bpsk/run_fig5.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import TABLE_I, Pe1_closed_form_eq39

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")


def main():
    gammas = np.linspace(1, 10, 20)
    fig, ax = plt.subplots(figsize=(7, 5))
    for regime, (a, b) in TABLE_I.items():
        for s, ls in [(0.1, '-'), (1.0, '--')]:
            vals = [Pe1_closed_form_eq39(g, a, b, s) for g in gammas]
            ax.semilogy(gammas, vals, ls, label=f'{regime}, σ={s}')
    ax.set_xlabel('|γ|'); ax.set_ylabel('Helstrom Error Probability (log)')
    ax.set_title('Hình 5 — 3 chế độ turbulence × 2 mức nhiễu (Table I)')
    ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "Fig5_reconstructed.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Da luu: {path}")


if __name__ == "__main__":
    main()
