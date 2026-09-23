""" generate_channel_samples.py -- Utility function to generate channel samples (Gamma-Gamma + Gaussian noise) shared across all experiments in experiments/bpsk/. """
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from model import sample_gamma_gamma

def generate_bpsk_channel_batch(alpha, beta, sigma, n_samples, seed=None, shared_noise=False):
    """Generate a batch of channel samples for the binary BPSK system |+gamma>, |-gamma>.
    shared_noise=True  : SHARED eps for both hypotheses (as used in paper eq.32)
    shared_noise=False : INDEPENDENT eps for both hypotheses (standard physical convention)
    Returns dict: {'I1':..., 'eps0':..., 'eps1':...}
    """
    rng = np.random.default_rng(seed)
    I1 = sample_gamma_gamma(alpha, beta, n_samples, rng)
    eps0 = rng.normal(0, sigma, n_samples)
    if shared_noise:
        eps1 = eps0.copy()
    else:
        eps1 = rng.normal(0, sigma, n_samples)
    return {"I1": I1, "eps0": eps0, "eps1": eps1}

def generate_bpsk_bitstream(alpha, beta, sigma, n_bits, seed=None):
    """Generate a random binary bitstream + corresponding channel for Figure 7-style simulations."""
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, n_bits)
    I1 = sample_gamma_gamma(alpha, beta, n_bits, rng)
    eps1 = rng.normal(0, sigma, n_bits)
    return bits, I1, eps1

if __name__ == "__main__":
    data = generate_bpsk_channel_batch(4.0, 1.9, 1.0, 5, seed=0)
    print("Channel data samples (5 samples):")
    for k, v in data.items():
        print(f"  {k}: {v}")
