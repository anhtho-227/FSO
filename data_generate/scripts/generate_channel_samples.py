"""
generate_channel_samples.py -- Hàm tiện ích sinh mẫu kênh (Gamma-Gamma + nhiễu Gauss)
dùng chung cho mọi thí nghiệm trong experiments/bpsk/.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from model import sample_gamma_gamma


def generate_bpsk_channel_batch(alpha, beta, sigma, n_samples, seed=None, shared_noise=False):
    """Sinh 1 lô mẫu kênh cho hệ BPSK nhị phân |+gamma>, |-gamma>.

    shared_noise=True  : eps DUNG CHUNG cho ca 2 gia thuyet (dung nhu eq.32 bai bao)
    shared_noise=False : eps DOC LAP cho 2 gia thuyet (quy uoc vat ly chuan)

    Tra ve dict: {'I1':..., 'eps0':..., 'eps1':...}
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
    """Sinh chuỗi bit nhị phân ngẫu nhiên + kênh tương ứng cho mô phỏng kiểu Hình 7."""
    rng = np.random.default_rng(seed)
    bits = rng.integers(0, 2, n_bits)
    I1 = sample_gamma_gamma(alpha, beta, n_bits, rng)
    eps1 = rng.normal(0, sigma, n_bits)
    return bits, I1, eps1


if __name__ == "__main__":
    data = generate_bpsk_channel_batch(4.0, 1.9, 1.0, 5, seed=0)
    print("Mau du lieu kenh (5 mau):")
    for k, v in data.items():
        print(f"  {k}: {v}")
