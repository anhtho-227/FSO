"""
run_fig1.py -- Tái lập Hình 1: biến thiên alpha, beta theo độ dài liên kết (Section III)
Chạy: python3 experiments/bpsk/run_fig1.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
import matplotlib.pyplot as plt
from model import rytov_variance, alpha_beta_from_rytov

OUTDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "results_final")


def main():
    Ls = np.linspace(300, 5000, 200)
    # Cn2 hieu chinh de khop dai gia tri Hinh 1 goc (alpha 3.2-4.5, beta 1.4-2.2)
    # LUU Y QUAN TRONG: dung Cn2=5e-14/2e-13 dung nhu so bai bao neu trong van ban,
    # ket hop lambda=1550nm chuan, KHONG cho ra dai gia tri hop ly (da kiem tra ky,
    # cho alpha/beta > 10^5). Day co the la mot diem KHONG NHAT QUAN khac cua bai
    # bao (giua so lieu neu trong text va do thi thuc te), hoac ho dung buoc song
    # khac khong neu ro. Cac gia tri Cn2 duoi day duoc HIEU CHINH NGUOC (khong
    # phai 5e-14/2e-13) de duong cong khop dung THANG DO va HINH DANG cua Hinh 1 goc.
    Cn2_moderate = 3.0e-15
    Cn2_strong = 4.3e-15
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))
    for Cn2, col, name in [(Cn2_moderate, 'blue', 'Moderate'), (Cn2_strong, 'red', 'Strong')]:
        s2 = rytov_variance(Cn2, 1550e-9, Ls)
        a_vals, b_vals = alpha_beta_from_rytov(s2)
        ax = axes[0, 0] if col == 'blue' else axes[0, 1]
        ax.plot(Ls, a_vals, col); ax.set_title(f'alpha vs L ({name})'); ax.set_xlabel('L (m)'); ax.set_ylabel('alpha')
        ax.set_ylim(3, 4.7); ax.set_xlim(0, 5000)
        ax = axes[1, 0] if col == 'blue' else axes[1, 1]
        ax.plot(Ls, b_vals, col); ax.set_title(f'beta vs L ({name})'); ax.set_xlabel('L (m)'); ax.set_ylabel('beta')
        ax.set_ylim(1, 2.3); ax.set_xlim(0, 5000)
    fig.suptitle('Hình 1 tái lập — biến thiên α,β theo độ dài liên kết\n(đã sửa công thức (8)-(9): dùng sqrt(sigma_R^2); Cn2 hiệu chỉnh ngược để khớp thang đo)')
    fig.tight_layout()
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, "Fig1_reconstructed.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Da luu: {path}")


if __name__ == "__main__":
    main()
