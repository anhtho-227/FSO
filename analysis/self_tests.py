""" 
self_tests.py 
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from model import (
    TABLE_I,
    helstrom_binary_pure,
    build_rho_matrix,
    helstrom_bound_as_written,
    helstrom_bound_corrected,
    rho_element_numeric_reference,
    rho_element_meijerG_corrected,
    Pe1_closed_form_eq39,
    Pe1_closed_form_robust,
)

def run_self_tests():
    print("=" * 88)
    print("SELF TESTS")
    print("=" * 88)
    
    alpha, beta = TABLE_I["moderate"]
    
    # --- Test 1: pure limit (sigma->0, I1=1 FIXED) matches Helstrom formula eq.5 ---
    print("\n[Test 1] Limit sigma->0, I1=1 FIXED (no turbulence): eq.(20) vs eq.(5)")
    from model.density_operator import rho_element_fixed_I
    sigma_small = 1e-3
    I1_fixed = 1.0
    N = 22
    for g in [0.5, 1.0, 1.5, 2.0]:
        rho0 = np.array([[rho_element_fixed_I(m, n, I1_fixed, +g, sigma_small) for n in range(N)] for m in range(N)])
        rho1 = np.array([[rho_element_fixed_I(m, n, I1_fixed, -g, sigma_small) for n in range(N)] for m in range(N)])
        pe = helstrom_bound_corrected(rho0, rho1)
        pe_pure = helstrom_binary_pure(g)
        ok = "OK" if abs(pe - pe_pure) < 1e-3 else "DEVIATED"
        print(f"  gamma={g} : model={pe:.6f} eq.(5)={pe_pure:.6f} [{ok}]")

    # --- Test 2: Meijer-G corrected matches numeric reference ---
    print("\n[Test 2] CASE 5: CORRECTED Meijer-G (k-term+1/2) vs numeric integration (on the fly)")
    sigma, g = 1.0, 2.0
    for (m, n) in [(0, 0), (1, 1), (2, 2), (0, 2)]:
        num = rho_element_numeric_reference(m, n, alpha, beta, +g, sigma)
        mei = rho_element_meijerG_corrected(m, n, alpha, beta, +g, sigma)
        dev = abs(num - mei) / abs(num) * 100 if num != 0 else float('nan')
        print(f"  (m,n)={(m,n)} : numeric={num:.8f} meijerG_corrected={mei:.8f} dev={dev:.4f} %")

    # --- Test 3: Meijer-G as-written for Tr != 1 and unreasonable Pe (bugs #2,#3) ---
    print("\n[Test 3] PAPER eq.(24) (bugs #2,#3) -> Tr(rho)!=1, Pe can be negative")
    N = 16
    for g in [1.0, 2.0]:
        rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="meijerG_as_written")
        rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="meijerG_as_written")
        pe = helstrom_bound_corrected(rho0, rho1)
        print(f"  gamma={g} : Tr(rho0)={np.trace(rho0):.4f} (should be 1)  Pe={pe:.4f} (can be negative/unreasonable)")

    # --- Test 4: eq.(26) without |.| always gives Pe=0.5 (bug #1) ---
    print("\n[Test 4] PAPER eq.(26) (without |.|) -> ALWAYS gives Pe=0.5")
    rho0 = build_rho_matrix(10, alpha, beta, +2.0, sigma, method="numeric")
    rho1 = build_rho_matrix(10, alpha, beta, -2.0, sigma, method="numeric")
    print(f"  eq.(26) without |.| -> Pe = {helstrom_bound_as_written(rho0, rho1):.6f} (always = 0.5, meaningless)")
    print(f"  eq.(25) with |.|    -> Pe = {helstrom_bound_corrected(rho0, rho1):.6f} (correct)")

    # --- Test 5: eq.(39) vs TRUE Helstrom -> deviates 3-8 times ---
    print("\n[Test 5] eq.(39) vs TRUE Helstrom limit (bug #4)")
    N = 20
    print(f" {'gamma':>6} | {'TRUE Pe (correct density-op)':>26} | {'Pe eq.(39)':>12} | {'ratio':>8} ")
    for g in [1.0, 2.0, 3.0]:
        rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="numeric")
        rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="numeric")
        pe_true = helstrom_bound_corrected(rho0, rho1)
        pe_39 = Pe1_closed_form_eq39(g, alpha, beta, sigma)
        print(f" {g:>6.1f} | {pe_true:>26.4f} | {pe_39:>12.4f} | {pe_true/pe_39:>7.2f} x")

    # --- Test 6: limit gamma->0 must be 0.5 ---
    print("\n[Test 6] Limit gamma->0 (Pe must = 0.5)")
    N = 18
    g = 1e-6
    rho0 = build_rho_matrix(N, alpha, beta, +g, sigma, method="numeric")
    rho1 = build_rho_matrix(N, alpha, beta, -g, sigma, method="numeric")
    print(f"  TRUE Pe (density-op)  = {helstrom_bound_corrected(rho0, rho1):.6f} (expected: 0.500000)")
    print(f"  Pe eq.(39) (robust)   = {Pe1_closed_form_robust(g, alpha, beta, sigma):.6f} (DOES NOT tend to 0.5)")

    print("\n" + "=" * 88)
    print("SELF-TEST SUITE COMPLETED")
    print("=" * 88)

if __name__ == "__main__":
    run_self_tests()
