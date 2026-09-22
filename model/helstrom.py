"""
helstrom.py -- LỚP 6: Helstrom detection qua trị riêng (Section VI-VII, eq.25-30)
LỖI #1 đã xác nhận bằng ảnh gốc trang 7: eq.(26)/(30) thiếu dấu trị tuyệt đối quanh Σλk
"""
import numpy as np


def helstrom_bound_as_written(rho0, rho1, p0=0.5, p1=0.5):
    """eq.(26)/(30) ĐÚNG Y VĂN BẢN: Pe = 1/2(1-Σλk) KHÔNG có |.|.
    Với p0=p1=0.5: Σλk = Tr(Δ) = p1-p0 = 0 LUÔN LUÔN -> Pe LUÔN = 0.5,
    vô giá trị đo lường (lỗi #1)."""
    Delta = p1*rho1 - p0*rho0
    Delta = (Delta+Delta.T)/2
    eig = np.linalg.eigvalsh(Delta)
    return 0.5*(1-np.sum(eig))


def helstrom_bound_corrected(rho0, rho1, p0=0.5, p1=0.5):
    """eq.(25) ĐÚNG: Pe = 1/2(1-Tr|Δ|) = 1/2(1-Σ|λk|)."""
    Delta = p1*rho1 - p0*rho0
    Delta = (Delta+Delta.T)/2
    eig = np.linalg.eigvalsh(Delta)
    return 0.5*(1-np.sum(np.abs(eig)))
