""" generate_table1.py -- Check/reproduce Table I (Section IX.C) from the Rytov model (Section III). 
Since Table I does not specify the exact C_n^2/wavelength/L used to derive the 3 rows, here we ONLY 
verify that the 3 pairs of (alpha, beta) in Table I lie on the same valid curve alpha(sigma_R^2), beta(sigma_R^2) 
(eq.8-9), by inversely solving for sigma_R^2 corresponding to each alpha and checking if the matching beta agrees. 
Run: python3 data_generate/scripts/generate_table1.py 
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
    print("Comparing Table I with the Rytov model (eq.8-9) -- NOTE: this formula")
    print("is non-monotonic (has a minimum around sigma_R^2 ~ 0.6-0.7), so a single alpha value")
    print("can correspond to two different sigma_R^2 values (weak / strong scintillation).\n")
    
    for regime, (a, b) in TABLE_I.items():
        s2 = invert_alpha(a)
        if s2 is not None:
            a_check, b_check = alpha_beta_from_rytov(s2)
            print(f"{regime:>8}: alpha={a}, beta={b}  ->  sigma_R^2~{s2:.4f} yields beta={b_check:.4f} "
                  f"(report table: {b})")
        else:
            print(f"{regime:>8}: could not find a matching sigma_R^2 for alpha={a} in the given range")

if __name__ == "__main__":
    main()
