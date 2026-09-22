"""
channel.py -- LỚP 2: Mô hình kênh Gamma-Gamma turbulence (Section III)  [ĐÃ XÁC NHẬN ĐÚNG]
"""
import numpy as np
from scipy.special import gamma as Gamma, kv as besselk


def rytov_variance(Cn2, wavelength, L, C=1.23):
    """sigma_R^2 = C.Cn2.k^(7/6).L^(11/6)   -- eq.(7)/(11)"""
    k = 2*np.pi/wavelength
    return C*Cn2*k**(7/6)*L**(11/6)


def alpha_beta_from_rytov(sigma2_R):
    """eq.(8)-(9) -- ĐÃ SỬA (xác nhận qua ảnh gốc trang 3-4, phóng to cực đại):
    công thức dùng sqrt(sigma_R^2), KHÔNG PHẢI sigma_R^2^(12/5) như công thức
    Andrews-Phillips "chuẩn" trong y văn thường thấy. Đây là điểm khác biệt riêng
    của bài báo này -- ĐÃ XÁC NHẬN BẰNG ẢNH, không phải suy đoán."""
    s = np.sqrt(sigma2_R)
    ra = 0.49*sigma2_R/(1+1.11*s)**(7/6)
    rb = 0.51*sigma2_R/(1+0.69*s)**(5/6)
    return 1/(np.exp(ra)-1), 1/(np.exp(rb)-1)


def gamma_gamma_pdf(i, alpha, beta):
    """f_I(i)  -- eq.(6)/(10)"""
    i = np.asarray(i, dtype=float)
    out = np.zeros_like(i)
    mask = i > 0
    pref = 2*(alpha*beta)**((alpha+beta)/2)/(Gamma(alpha)*Gamma(beta))
    out[mask] = pref*i[mask]**((alpha+beta)/2-1)*besselk(alpha-beta, 2*np.sqrt(alpha*beta*i[mask]))
    return out


def sample_gamma_gamma(alpha, beta, n, rng=None):
    """Sinh n mẫu I ~ Gamma-Gamma(alpha,beta), E[I]=1."""
    rng = rng or np.random.default_rng()
    X = rng.gamma(shape=alpha, scale=1.0/alpha, size=n)
    Y = rng.gamma(shape=beta, scale=1.0/beta, size=n)
    return X*Y
