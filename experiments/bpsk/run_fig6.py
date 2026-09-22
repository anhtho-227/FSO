"""
run_fig6.py -- Tái lập Hình 6: Pe, Pguess ideal vs thực tế (Section IX.D, eq.41-43)
Chạy: python3 experiments/bpsk/run_fig6.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import Pe_ideal, Pguess_ideal, Pe_real, Pguess_real

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")


def main():
    alpha, beta, sigma = 4.2, 1.4, 0.1  # strong turbulence, xem trang 10-11
    gammas = np.linspace(0.05, 3.0, 25)
    pe_i = [Pe_ideal(g) for g in gammas]
    pe_r = [Pe_real(g, alpha, beta, sigma) for g in gammas]
    pg_ideal = [1-p for p in pe_i]
    pg_real = [1-p for p in pe_r]

    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.plot(gammas, pg_ideal, 'b-^', label='Pguess (Ideal)', markersize=4)
    ax.plot(gammas, pg_real, 'r-v', label='Pguess (Turbulence+Noise)', markersize=4)
    ax.plot(gammas, pe_i, 'b--o', label='Pe (Ideal)', markersize=3, markerfacecolor='none')
    ax.plot(gammas, pe_r, 'r--s', label='Pe (Turbulence+Noise)', markersize=3, markerfacecolor='none')
    for thr in [0.75, 0.80, 0.90]:
        ax.axhline(thr, color='k', linestyle=':', linewidth=0.8)
    ax.set_xlabel('γ (Field Amplitude)'); ax.set_ylabel('Probability')
    ax.set_title('Hình 6 — Relay detection & Eavesdropper confidence')
    ax.grid(True, alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "Fig6_reconstructed.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Da luu: {path}")


if __name__ == "__main__":
    main()
