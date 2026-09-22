"""
density_operator.py
LỚP 4: Toán tử mật độ, I1 CỐ ĐỊNH (Section V.A-C, eq.13-20)      [ĐÃ XÁC NHẬN ĐÚNG]
LỚP 5: Toán tử mật độ, TRUNG BÌNH Gamma-Gamma (Appendix A-B, eq.21-24/60)
       eq.(22)-(24) CÓ 2 LỖI đã xác nhận bằng ảnh gốc PDF:
         #2 thiếu điều kiện "k chẵn" (trang 13, eq.57 có / eq.22 không)
         #3 thiếu hệ số ½ kế thừa từ eq.(60) (trang 13/7)
"""
import math
import numpy as np
import mpmath as mp
from scipy.special import gamma as Gamma, comb
from scipy import integrate

from .channel import gamma_gamma_pdf

mp.mp.dps = 25


# ---------------------------------------------------------------------------
# LỚP 4
# ---------------------------------------------------------------------------

def rho_element_fixed_I(m, n, I1, psi_i, sigma):
    """[rho_i(I1)]_{m,n}  -- eq.(20)/(55).
    Ảnh gốc trang 6,12 XÁC NHẬN có điều kiện "k chẵn" -- đã kiểm chứng khớp
    tuyệt đối với tích phân trực tiếp của eq.(19)."""
    A = (1+2*sigma**2)/(2*sigma**2)
    pref = 1/(np.sqrt(2*np.pi)*sigma)*np.exp(-I1*psi_i**2/(1+2*sigma**2))
    total = 0.0
    for k in range(0, m+n+1):
        if k % 2 != 0:
            continue
        total += comb(m+n, k)*Gamma((k+1)/2)/(A**((k+1)/2)) \
                 * (np.sqrt(I1)*psi_i/(1+2*sigma**2))**(m+n-k)
    return pref*total/math.sqrt(math.factorial(m)*math.factorial(n))


# ---------------------------------------------------------------------------
# LỚP 5
# ---------------------------------------------------------------------------

def rho_element_numeric_reference(m, n, alpha, beta, psi_i, sigma):
    """'TRỌNG TÀI' -- tích phân số trực tiếp eq.(21): ∫ eq.(20).f_Gamma-Gamma di.
    Đáng tin cậy nhất (không qua Meijer-G, tránh lỗi tham số hàm đặc biệt).
    Đã kiểm chứng: giới hạn sigma->0 khớp tuyệt đối eq.(5); khớp 0.0000% với
    Meijer-G đã sửa (xem analysis/self_tests.py, Test 2)."""
    f = lambda i: rho_element_fixed_I(m, n, i, psi_i, sigma) * gamma_gamma_pdf(np.array([i]), alpha, beta)[0]
    val, _ = integrate.quad(f, 1e-9, 100, limit=200)
    return val


def rho_element_meijerG(m, n, alpha, beta, psi_i, sigma, k_even=True, half_factor=True):
    """[rho_i]_{m,n} qua Meijer-G -- eq.(22)-(24)/(60).

    Quy ước tham số Meijer-G ĐÚNG: G^{2,1}_{1,2}(z | a1 ; b1,b2)
        -> mp.meijerg([[a1],[]], [[b1,b2],[]], z)   (n=1,p=1 ; m=2,q=2)

    k_even=False, half_factor=False -> ĐÚNG Y VĂN BẢN eq.(24) (có lỗi #2+#3)
    k_even=True,  half_factor=True  -> ĐÃ SỬA (khớp 0.0000% với trọng tài)
    """
    m_, n_, alpha_, beta_, psi_, sigma_ = [mp.mpf(x) for x in (m, n, alpha, beta, psi_i, sigma)]
    pref = (2/(sigma_*mp.sqrt(2*mp.pi))) * (alpha_*beta_)**((alpha_+beta_)/2) / (mp.gamma(alpha_)*mp.gamma(beta_))
    total = mp.mpf(0)
    for k in range(m+n+1):
        if k_even and k % 2 != 0:
            continue
        kf = mp.mpf(k)
        term = mp.binomial(m_+n_, kf)/mp.sqrt(mp.factorial(m_)*mp.factorial(n_))
        term *= (1/(1+2*sigma_**2))**((m_+n_-kf-alpha_-beta_)/2)
        term *= (2*sigma_**2/(1+2*sigma_**2))**((kf+1)/2)
        term *= mp.gamma((kf+1)/2)
        term *= abs(psi_)**(-m_-n_-alpha_-beta_+kf) * psi_**(m_+n_-kf)
        z = alpha_*beta_*(1+2*sigma_**2)/psi_**2
        a1 = 1-(m_+n_-kf+alpha_+beta_)/2
        b1, b2 = (alpha_-beta_)/2, -(alpha_-beta_)/2
        Gval = mp.meijerg([[a1], []], [[b1, b2], []], z)
        if half_factor:
            Gval *= mp.mpf('0.5')
        total += term*Gval
    return float(pref*total)


def rho_element_meijerG_as_written(m, n, alpha, beta, psi_i, sigma):
    """eq.(24) ĐÚNG Y VĂN BẢN GỐC (có lỗi #2 và #3) -- Tr(rho) != 1, có thể cho
    Helstrom bound ÂM (vô lý) nếu dùng để tính hệ thống."""
    return rho_element_meijerG(m, n, alpha, beta, psi_i, sigma, k_even=False, half_factor=False)


def rho_element_meijerG_corrected(m, n, alpha, beta, psi_i, sigma):
    """eq.(24) ĐÃ SỬA lỗi #2 (k chẵn) + lỗi #3 (hệ số ½)."""
    return rho_element_meijerG(m, n, alpha, beta, psi_i, sigma, k_even=True, half_factor=True)


def build_rho_matrix(N, alpha, beta, psi_i, sigma, method="numeric"):
    """Dựng ma trận mật độ NxN cắt cụt Fock basis.
    method: 'numeric' (trọng tài, khuyến nghị) | 'meijerG_corrected' | 'meijerG_as_written'
    """
    fn = {
        "numeric": rho_element_numeric_reference,
        "meijerG_corrected": rho_element_meijerG_corrected,
        "meijerG_as_written": rho_element_meijerG_as_written,
    }[method]
    rho = np.zeros((N, N))
    for m in range(N):
        for n in range(m, N):
            v = fn(m, n, alpha, beta, psi_i, sigma)
            rho[m, n] = v
            rho[n, m] = v
    return rho
