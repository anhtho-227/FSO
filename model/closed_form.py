"""
closed_form.py -- LỚP 7: Công thức đóng xấp xỉ (Section VIII, Appendix C-E, eq.31-40)

LỖI #4 (nghiêm trọng nhất, đã xác nhận bằng tính toán độc lập nhiều lớp):
quy ước "nhiễu chia sẻ, đảo dấu" ở eq.(32) khiến TOÀN BỘ nhánh này đánh giá
THẤP giới hạn Helstrom thật 3-8 lần. Bản thân đại số từ (32)->(39) là ĐÚNG,
vấn đề nằm ở điểm xuất phát vật lý (eq.32), không phải các bước suy diễn.
"""
import numpy as np
import mpmath as mp
from scipy.special import gamma as Gamma, hyp1f1
from scipy import integrate

from .channel import gamma_gamma_pdf

mp.mp.dps = 25


def overlap_shared_noise_as_written(I1, gamma_, eps1):
    """<psi0|psi1> = exp(-2|sqrt(I1)gamma+eps1|^2)  -- eq.(34), ĐÚNG Y VĂN BẢN
    (dùng psi1 = -(sqrt(I1)gamma+eps1), nhiễu "chia sẻ + đảo dấu" -- gốc rễ lỗi #4)."""
    a = np.sqrt(I1)*gamma_+eps1
    return np.exp(-2*np.abs(a)**2)


def Pe_exact_conditional_eq36(I1, gamma_, eps1):
    """eq.(36): Pe = 1/2(1-sqrt(1-|overlap|^2)), CHƯA xấp xỉ Taylor."""
    c2 = np.abs(overlap_shared_noise_as_written(I1, gamma_, eps1))**2
    return 0.5*(1-np.sqrt(np.maximum(0, 1-c2)))


def Pe_approx_conditional_eq37(I1, gamma_, eps1):
    """eq.(37): Pe ≈ 1/4.exp(-4|sqrt(I1)gamma+eps1|^2)
    (xấp xỉ (1-x)^0.5≈1-x/2, chỉ hợp lệ khi 4|.|^2>>1)."""
    a = np.sqrt(I1)*gamma_+eps1
    return 0.25*np.exp(-4*a**2)


def Pe1_averaged_over_noise_eq38(I1, gamma_, sigma):
    """eq.(38) -- công thức đóng của E_eps1[eq.37]. ĐÃ KIỂM CHỨNG ĐÚNG."""
    return np.exp(-4*I1*gamma_**2/(1+8*sigma**2)) / (4*np.sqrt(1+8*sigma**2))


def Pe1_closed_form_eq39(gamma_, alpha, beta, sigma):
    """eq.(39)/(76) -- ĐÚNG THEO ẢNH GỐC. Về ĐẠI SỐ, suy diễn đúng từ eq.(38);
    nhưng GIÁ TRỊ thấp hơn giới hạn Helstrom THẬT khoảng 3-8 lần (lỗi #4)."""
    b = 4*gamma_**2/(1+8*sigma**2)
    pref = 0.25*np.sqrt(1/(1+8*sigma**2))*(alpha*beta)**beta/(Gamma(alpha)*Gamma(beta))
    arg = alpha*beta/b
    term1 = (b**(alpha-beta))*Gamma(alpha-beta)*Gamma(beta) * hyp1f1(beta, 1-(alpha-beta), arg)
    term2 = (alpha*beta)**(alpha-beta)*Gamma(beta-alpha)*Gamma(alpha) * hyp1f1(alpha, 1+(alpha-beta), arg)
    return pref*b**(-alpha)*(term1+term2)


def Pe1_closed_form_hp(gamma_, alpha, beta, sigma):
    """Bản mpmath độ chính xác cao của eq.(39) -- dùng khi gamma rất nhỏ làm
    scipy.hyp1f1 tràn số (huỷ triệt tiêu giữa 2 số hạng lớn)."""
    g_, a_, b_, s_ = [mp.mpf(x) for x in (gamma_, alpha, beta, sigma)]
    bb = 4*g_**2/(1+8*s_**2)
    pref = mp.mpf('0.25')*mp.sqrt(1/(1+8*s_**2))*(a_*b_)**b_/(mp.gamma(a_)*mp.gamma(b_))
    arg = a_*b_/bb
    term1 = (bb**(a_-b_))*mp.gamma(a_-b_)*mp.gamma(b_) * mp.hyp1f1(b_, 1-(a_-b_), arg)
    term2 = (a_*b_)**(a_-b_)*mp.gamma(b_-a_)*mp.gamma(a_) * mp.hyp1f1(a_, 1+(a_-b_), arg)
    return float(pref*bb**(-a_)*(term1+term2))


def Pe1_closed_form_robust(gamma_, alpha, beta, sigma):
    """Bản ỔN ĐỊNH SỐ HỌC của eq.(39): tích phân số trực tiếp theo I thay vì
    tổng 2 số hạng hypergeometric (vốn huỷ triệt tiêu số học ở gamma nhỏ)."""
    def L_of_I(i):
        return np.exp(-4*i*gamma_**2/(1+8*sigma**2)) / (4*np.sqrt(1+8*sigma**2))
    f = lambda i: L_of_I(i)*gamma_gamma_pdf(np.array([i]), alpha, beta)[0]
    val, _ = integrate.quad(f, 1e-8, 300, limit=400)
    return min(val, 1.0)
