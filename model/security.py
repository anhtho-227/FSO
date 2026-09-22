"""
security.py -- LỚP 8: Bảo mật, xác suất đoán đúng của relay (Section IX.D, eq.41-43)
"""
from .quantum_states import helstrom_binary_pure
from .closed_form import Pe1_closed_form_robust


def Pe_ideal(gamma_):
    """eq.(41)"""
    return helstrom_binary_pure(gamma_)


def Pguess_ideal(gamma_):
    """eq.(42)"""
    return 1-Pe_ideal(gamma_)


def Pe_real(gamma_, alpha, beta, sigma):
    """Dùng eq.(37)/(39) closed-form, bản ổn định số học -- theo mô tả trang 10-11."""
    return Pe1_closed_form_robust(gamma_, alpha, beta, sigma)


def Pguess_real(gamma_, alpha, beta, sigma):
    """eq.(43)"""
    return 1-Pe_real(gamma_, alpha, beta, sigma)
