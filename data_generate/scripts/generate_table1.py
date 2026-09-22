"""
generate_table1.py -- Kiểm tra/tái tạo Table I (Section IX.C) từ mô hình Rytov (Section III).

Table I không cho biết C_n^2/wavelength/L cụ thể dùng để suy ra 3 hàng, nên ở đây
ta CHỈ kiểm tra rằng bộ 3 cặp (alpha,beta) trong Table I nằm trên cùng một đường
cong alpha(sigma_R^2), beta(sigma_R^2) hợp lệ (eq.8-9), bằng cách dò ngược sigma_R^2
ứng với từng alpha rồi kiểm tra beta tương ứng có khớp không.

Chạy: python3 data_generate/scripts/generate_table1.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

import numpy as np
from scipy.optimize import brentq
from model import TABLE_I
from model.channel import alpha_beta_from_rytov


def invert_alpha(alpha_target, lo=1e-6, hi=100):
    f = lambda s2: alpha_beta_from_rytov(s2)[0] - alpha_target
    try:
        return brentq(f, lo, hi)
    except ValueError:
        return None


def main():
    print("Doi chieu Table I voi mo hinh Rytov (eq.8-9) -- LUU Y: cong thuc nay")
    print("khong don dieu (co cuc tieu quanh sigma_R^2~0.6-0.7), nen 1 gia tri alpha")
    print("co the ung voi 2 gia tri sigma_R^2 khac nhau (nhanh yeu / nhanh manh).\n")
    for regime, (a, b) in TABLE_I.items():
        s2 = invert_alpha(a)
        if s2 is not None:
            a_check, b_check = alpha_beta_from_rytov(s2)
            print(f"{regime:>8}: alpha={a}, beta={b}  ->  sigma_R^2~{s2:.4f} suy ra beta={b_check:.4f} "
                  f"(bang bao cao: {b})")
        else:
            print(f"{regime:>8}: khong tim duoc sigma_R^2 khop alpha={a} trong khoang do")


if __name__ == "__main__":
    main()
