"""
run_fig4.py -- Tái lập Hình 4: Direct vs Relay, weak turbulence (Section IX.B)
Chạy: python3 experiments/bpsk/run_fig4.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import DIRECT_LINKS, RELAY_LINKS, Pe1_closed_form_eq39, end_to_end_error

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")


def main():
    sigma = 0.1
    gammas = np.linspace(1, 10, 20)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    colors = plt.cm.tab10(np.linspace(0, 1, 4))
    for idx, L in enumerate(DIRECT_LINKS.keys()):
        a_d, b_d = DIRECT_LINKS[L]
        a_r, b_r = RELAY_LINKS[L // 2]
        pe_direct = [Pe1_closed_form_eq39(g, a_d, b_d, sigma) for g in gammas]
        pe_hop = [Pe1_closed_form_eq39(g, a_r, b_r, sigma) for g in gammas]
        pe_relay = [end_to_end_error(p, p) for p in pe_hop]
        ax.semilogy(gammas, pe_direct, '*-', color=colors[idx], label=f'Direct, L={L}m')
        ax.semilogy(gammas, pe_relay, 'o--', color=colors[idx], markerfacecolor='none', label=f'Relay, L={L//2}m x2')
    ax.set_xlabel('|γ|'); ax.set_ylabel('Helstrom Error Probability (log)')
    ax.set_title('Hình 4 — Direct vs Relay, weak turbulence, σ=0.1\n(dùng eq.39 — lưu ý: thấp hơn giá trị thật ~3-8 lần, xem lỗi #4)')
    ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(fontsize=7)
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "Fig4_reconstructed.png")
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"Da luu: {path}")


if __name__ == "__main__":
    main()
