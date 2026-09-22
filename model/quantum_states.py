"""
quantum_states.py -- LỚP 1: Trạng thái lượng tử cơ bản (Section II)  [ĐÃ XÁC NHẬN ĐÚNG]
"""
import numpy as np


def coherent_overlap(a, b):
    """<a|b> = exp(-|a|^2/2 - |b|^2/2 + a*.b)   -- eq.(2)"""
    return np.exp(-0.5*abs(a)**2 - 0.5*abs(b)**2 + np.conj(a)*b)


def helstrom_binary_pure(alpha_amp):
    """P_e^min = 1/2(1-sqrt(1-exp(-4|alpha|^2)))  cho |+alpha> vs |-alpha>  -- eq.(5)"""
    return 0.5*(1-np.sqrt(1-np.exp(-4*abs(alpha_amp)**2)))
