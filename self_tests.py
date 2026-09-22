"""
self_tests.py -- LỚP 9: BỘ TỰ KIỂM CHỨNG

Chạy: python3 -m analysis.self_tests   (từ thư mục gốc quantum_relay_reconstruction/)
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from model import (
    TABLE_I, helstrom_binary_pure, build_rho_matrix,
    helstrom_bound_as_written, helstrom_bound_corrected,
    rho_element_numeric_reference, rho_element_meijerG_corrected,
    Pe1_closed_form_eq39, Pe1_closed_form_robust,
)


def run_self_tests():
    print("="*88)
    print("BO TU KIEM CHUNG -- quantum_relay_reconstruction")
    print("="*88)
    alpha, beta = TABLE_I["moderate"]

    # --- Test 1: gioi han thuan khiet (sigma->0, I1=1 CO DINH) khop cong thuc Helstrom eq.5 ---
    print("\n[Test 1] Gioi han sigma->0, I1=1 CO DINH (khong turbulence): eq.(20) vs eq.(5)")
    from model.density_operator import rho_element_fixed_I
    sigma_small = 1e-3
    I1_fixed = 1.0
    N = 22
    for g in [0.5, 1.0, 1.5, 2.0]:
        rho0 = np.array([[rho_element_fixed_I(m, n, I1_fixed, +g, sigma_small) for n in range(N)] for m in range(N)])
        rho1 = np.array([[rho_element_fixed_I(m, n, I1_fixed, -g, sigma_small) for n in range(N)] for m in range(N)])
        pe = helstrom_bound_corrected(rho0, rho1)
        pe_pure = helstrom_binary_pure(g)
        ok = "OK" if abs(pe-pe_pure) < 1e-3 else "LECH"
        print(f"  gamma={g}: model={pe:.6f}  eq.(5)={pe_pure:.6f}  [{ok}]")

    # --- Test 2: Meijer-G corrected khop numeric reference ---
    print("\n[Test 2] LOP 5: Meijer-G DA SUA (k-chan+1/2) vs tich phan so (trong tai)")
    sigma, g = 1.0, 2.0
    for (m, n) in [(0, 0), (1, 1), (2, 2), (0, 2)]:
        num = rho_element_numeric_reference(m, n, alpha, beta, +g, sigma)
        mei = rho_element_meijerG_corrected(m, n, alpha, beta, +g, sigma)
        lech = abs(num-mei)/abs(num)*100 if num != 0 else float('nan')
        print(f"  (m,n)={(m,n)}: numeric={num:.8f}  meijerG_corrected={mei:.8f}  lech={lech:.4f}%")

    # --- Test 3: Meijer-G as-written cho Tr != 1 va Pe vo ly (loi #2,#3) ---
    print("\n[Test 3] eq.(24) DUNG Y VAN BAN (loi #2,#3) -> Tr(rho)!=1, Pe co the am")
    N = 16
    for g in [1.0, 2.0]:
        rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="meijerG_as_written")
        rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="meijerG_as_written")
        pe = helstrom_bound_corrected(rho0, rho1)
        print(f"  gamma={g}: Tr(rho0)={np.trace(rho0):.4f} (nen =1)  Pe={pe:.4f} (co the am/vo ly)")

    # --- Test 4: eq.(26) khong |.| luon cho Pe=0.5 (loi #1) ---
    print("\n[Test 4] eq.(26) DUNG Y VAN BAN (khong |.|) -> LUON cho Pe=0.5")
    rho0 = build_rho_matrix(10, alpha, beta, +2.0, sigma, method="numeric")
    rho1 = build_rho_matrix(10, alpha, beta, -2.0, sigma, method="numeric")
    print(f"  eq.(26) khong |.|  -> Pe = {helstrom_bound_as_written(rho0, rho1):.6f}  (luon = 0.5, vo nghia)")
    print(f"  eq.(25) co |.|     -> Pe = {helstrom_bound_corrected(rho0, rho1):.6f}  (dung)")

    # --- Test 5: eq.(39) vs Helstrom THAT -> lech 3-8 lan ---
    print("\n[Test 5] eq.(39) vs gioi han Helstrom THAT (loi #4)")
    N = 20
    print(f"  {'gamma':>6} | {'Pe THAT (density-op dung)':>26} | {'Pe eq.(39)':>12} | {'ty le':>8}")
    for g in [1.0, 2.0, 3.0]:
        rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="numeric")
        rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="numeric")
        pe_true = helstrom_bound_corrected(rho0, rho1)
        pe_39 = Pe1_closed_form_eq39(g, alpha, beta, sigma)
        print(f"  {g:>6.1f} | {pe_true:>26.4f} | {pe_39:>12.4f} | {pe_true/pe_39:>7.2f}x")

    # --- Test 6: gioi han gamma->0 phai la 0.5 ---
    print("\n[Test 6] Gioi han gamma->0 (Pe phai = 0.5)")
    N = 18
    g = 1e-6
    rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="numeric")
    rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="numeric")
    print(f"  Pe THAT (density-op)  = {helstrom_bound_corrected(rho0, rho1):.6f}  (ky vong: 0.500000)")
    print(f"  Pe eq.(39) (on dinh)  = {Pe1_closed_form_robust(g, alpha, beta, sigma):.6f}  (KHONG tien ve 0.5)")

    print("\n" + "="*88)
    print("HOAN TAT BO TU KIEM CHUNG")
    print("="*88)


if __name__ == "__main__":
    run_self_tests()
