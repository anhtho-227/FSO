import numpy as np
from scipy.special import gamma as gamma_func, hyp1f1, comb
import mpmath
import sys
sys.stdout.reconfigure(encoding='utf-8')
mpmath.mp.dps = 25

def helstrom_closed_form_strict(gamma_amp, alpha, beta, sigma):
    """Tính chính xác theo Phương trình 39, không sửa xấp xỉ Taylor"""
    term_gamma2 = (4 * gamma_amp**2) / (1 + 8 * sigma**2)
    ab = alpha * beta
    coeff = 0.25 * np.sqrt(1 / (1 + 8 * sigma**2)) * (ab**beta) / \
            (gamma_func(alpha) * gamma_func(beta)) * (term_gamma2**(-alpha))
    part1 = (term_gamma2**(alpha - beta)) * gamma_func(alpha - beta) * gamma_func(beta) * \
            hyp1f1(beta, 1 - (alpha - beta), ab / term_gamma2)
    part2 = (ab**(alpha - beta)) * gamma_func(beta - alpha) * gamma_func(alpha) * \
            hyp1f1(alpha, 1 + (alpha - beta), ab / term_gamma2)
    return coeff * (part1 + part2)

def gamma_gamma_pdf_params(alpha, beta):
    return 2 * (alpha * beta)**((alpha + beta) / 2) / (gamma_func(alpha) * gamma_func(beta))

def calc_rho_element_strict(m, n, gamma_amp, alpha, beta, sigma, state_sign):
    """Tính Phương trình 24 với mọi k, bỏ qua các điều kiện chẵn/lẻ"""
    psi = state_sign * gamma_amp
    term1 = 1 / (sigma * np.sqrt(2 * np.pi))
    term2 = gamma_gamma_pdf_params(alpha, beta)
    sigma_factor = 1 + 2 * sigma**2
    sum_val = mpmath.mpf(0)
    
    for k in range(m + n + 1):
        comb_coeff = comb(m + n, k) / np.sqrt(float(mpmath.fac(m) * mpmath.fac(n)))
        pow1 = (1 / sigma_factor)**((m + n - k - alpha - beta) / 2)
        pow2 = ((2 * sigma**2) / sigma_factor)**((k + 1) / 2)
        gamma_k = float(mpmath.gamma((k + 1) / 2)) # Lỗi toán học ở đây khi k lẻ
        psi_pow = (abs(psi)**(-m - n - alpha - beta + k)) * (psi**(m + n - k))
        
        z = (alpha * beta * sigma_factor) / (psi**2)
        a1 = 1 - (m + n - k + alpha + beta) / 2
        b1, b2 = (alpha - beta) / 2, -(alpha - beta) / 2
        
        meijer_g = mpmath.meijerg([[a1], []], [[b1, b2], []], z)
        sum_val += comb_coeff * pow1 * pow2 * gamma_k * psi_pow * meijer_g
        
    return float(term1 * term2 * sum_val)

def verify_paper_flaws():
    
    alpha, beta, sigma = 4.0, 1.9, 1.0
    gamma_zero = 0.0001 
    gamma_test = 2.0
    N = 4 
    
    pe_zero = helstrom_closed_form_strict(gamma_zero, alpha, beta, sigma)
    print(f"[TEST 1 - Phương trình 39] Giới hạn tín hiệu tiến về 0:")
    print(f"  -> Trạng thái thực tế: 2 tín hiệu trùng nhau (Lỗi 50%)")
    print(f"  -> Kết quả từ công thức bài báo: {pe_zero:.4f} ({(pe_zero*100):.2f}%)")
    
    rho_0 = np.array([[calc_rho_element_strict(m, n, gamma_test, alpha, beta, sigma, 1) for n in range(N)] for m in range(N)])
    trace_rho_0 = np.trace(rho_0)
    print(f"[TEST 2 - Phương trình 24] Tổng vết ma trận mật độ (Trace):")
    print(f"  -> Trạng thái thực tế: Tổng xác suất trên toàn không gian phải bằng 1.0")
    print(f"  -> Kết quả từ công thức bài báo: {trace_rho_0:.6f}")
    
    rho_1 = np.array([[calc_rho_element_strict(m, n, gamma_test, alpha, beta, sigma, -1) for n in range(N)] for m in range(N)])
    Delta = 0.5 * rho_1 - 0.5 * rho_0
    eigenvalues = np.linalg.eigvalsh(Delta)

if __name__ == "__main__":
    verify_paper_flaws()