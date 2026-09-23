""" run_fig3.py -- Reconstruct Figure 3 (Section IX.A), 2 modes:
--mode match     : matches EXACT MAGNITUDE as the original paper (uses "sample-by-sample conditioning, independent noise" calculation -- quantity B -- closely matches eq.39)
--mode rigorous  : TRUE Helstrom (density-operator with bugs #2+#3 fixed, standard trace-norm) -- shows a 3-8x deviation compared to eq.(39) (bug #4)
"""
import sys, os, argparse
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import TABLE_I, Pe1_closed_form_eq39, Pe_exact_conditional_eq36, build_rho_matrix, helstrom_bound_corrected
from data_generate.scripts.generate_channel_samples import generate_bpsk_channel_batch

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")

def fock_truncated_pure_state_pe(gamma_, alpha, beta, sigma, N, n_samples, seed):
    """Calculation method for 'quantity B': sample-by-sample conditioning (I1, INDEPENDENT noise for two hypotheses),
    constructing an N-dimensional Fock-truncated pure state, calculating Helstrom for that pair of pure states, and averaging --
    this is the method that makes the magnitude match Figure 3 of the paper (see ERRORS.md, section "Why V2 didn't detect bug #4")."""
    data = generate_bpsk_channel_batch(alpha, beta, sigma, n_samples, seed=seed, shared_noise=False)
    I1, eps_a, eps_b = data["I1"], data["eps0"], data["eps1"]
    beta0 = np.sqrt(I1)*gamma_ + eps_a
    beta1 = -np.sqrt(I1)*gamma_ + eps_b
    
    n_idx = np.arange(N)
    log_fact = np.concatenate([[0.0], np.cumsum(np.log(np.arange(1, N)))]) if N > 1 else np.array([0.0])
    
    def fock_vec(b):
        sign = np.sign(b)[:, None] ** n_idx[None, :]
        log_amp = -0.5*b[:, None]**2 + n_idx[None, :]*np.log(np.abs(b)[:, None]+1e-300) - 0.5*log_fact[None, :]
        return sign * np.exp(log_amp)
        
    v0 = fock_vec(beta0); v0 = v0/np.linalg.norm(v0, axis=1, keepdims=True)
    v1 = fock_vec(beta1); v1 = v1/np.linalg.norm(v1, axis=1, keepdims=True)
    overlap = np.sum(v0*v1, axis=1)
    pe = 0.5*(1-np.sqrt(np.maximum(0, 1-overlap**2)))
    return pe.mean()

def run_match(alpha, beta, sigma, gammas, outpath):
    N_list = [4, 8, 12, 16, 20, 24, 28, 32]
    pe_closed = np.array([Pe1_closed_form_eq39(g, alpha, beta, sigma) for g in gammas])
    data = generate_bpsk_channel_batch(alpha, beta, sigma, 400_000, seed=2, shared_noise=True)
    pe_ensemble = np.array([Pe_exact_conditional_eq36(data["I1"], g, data["eps0"]).mean() for g in gammas])
    
    pe_by_N = {}
    for N in N_list:
        pe_by_N[N] = np.array([fock_truncated_pure_state_pe(g, alpha, beta, sigma, N, 150_000, 100+N) for g in gammas])
        
    fig, ax = plt.subplots(figsize=(7, 5.2))
    for N, pe_N in pe_by_N.items():
        ax.semilogy(gammas, pe_N, '--', color='tab:blue', alpha=0.6, linewidth=1)
        ax.annotate(f'{N} Fock basis', xy=(gammas[1], pe_N[1]*1.15), fontsize=6.5, color='tab:blue')
        
    ax.semilogy(gammas, pe_closed, '-o', color='tab:red', label='Closed-form eq.(39)', markerfacecolor='none')
    ax.semilogy(gammas, pe_ensemble, '-^', color='tab:cyan', label='Ensemble MC eq.(36)')
    ax.set_xlabel('|γ|'); ax.set_ylabel('Helstrom Error Probability (log)')
    ax.set_title(f'Figure 3 [match] — α={alpha}, β={beta}, σ={sigma}')
    ax.set_ylim(1e-4, 1e0); ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outpath, dpi=160); plt.close(fig)

def run_rigorous(alpha, beta, sigma, gammas, outpath, N=24):
    pe_closed = np.array([Pe1_closed_form_eq39(g, alpha, beta, sigma) for g in gammas])
    pe_true = []
    for g in gammas:
        rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="numeric")
        rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="numeric")
        pe_true.append(helstrom_bound_corrected(rho0, rho1))
    pe_true = np.array(pe_true)
    
    fig, ax = plt.subplots(figsize=(7, 5.2))
    ax.semilogy(gammas, pe_closed, '-o', color='tab:red', label='Closed-form eq.(39) [Section VIII]')
    ax.semilogy(gammas, pe_true, '-s', color='black', linewidth=2, label=f'TRUE Helstrom (fixed density-op, N={N})')
    ax.set_xlabel('|γ|'); ax.set_ylabel('Helstrom Error Probability (log)')
    ax.set_title(f'Figure 3 [rigorous] — α={alpha}, β={beta}, σ={sigma}\n(shows bug #4: eq.39 is 3-8x lower)')
    ax.set_ylim(1e-4, 1e0); ax.grid(True, which='both', ls=':', alpha=0.5)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(outpath, dpi=160); plt.close(fig)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["match", "rigorous"], default="match")
    args = parser.parse_args()
    
    alpha, beta = TABLE_I["moderate"]
    sigma = 1.0
    gammas = np.array([1.0, 1.3, 1.7, 2.1, 2.4, 2.8, 3.2, 3.6, 4.0])
    
    os.makedirs(OUTDIR, exist_ok=True)
    if args.mode == "match":
        outpath = os.path.join(OUTDIR, "Fig3_match_paper.png")
        run_match(alpha, beta, sigma, gammas, outpath)
    else:
        outpath = os.path.join(OUTDIR, "Fig3_rigorous.png")
        run_rigorous(alpha, beta, sigma, gammas, outpath)
        
    print(f"Saved: {outpath}")

if __name__ == "__main__":
    main()
