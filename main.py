import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as gamma_func, hyp1f1, comb
import mpmath
import sys
sys.stdout.reconfigure(encoding='utf-8')
mpmath.mp.dps = 25

# =====================================================================
# KHỐI 1: MÔ HÌNH KÊNH VẬT LÝ VÀ NHIỄU LOẠN KHÍ QUYỂN
# =====================================================================
def calculate_gamma_gamma_params(L, Cn2, wavelength_nm=1550):
    """Tính thông số alpha, beta từ khoảng cách L qua phương sai Rytov (Phương trình 7, 8, 9)."""
    k = 2 * np.pi / (wavelength_nm * 1e-9)
    sigma_R2 = 1.23 * Cn2 * (k**(7/6)) * (L**(11/6)) 
    
    alpha = (np.exp(0.49 * sigma_R2 / (1 + 1.11 * np.sqrt(sigma_R2))**(1.2)) - 1)**(-1)
    beta = (np.exp(0.51 * sigma_R2 / (1 + 0.69 * np.sqrt(sigma_R2))**(1.2)) - 1)**(-1) 
    return alpha, beta, sigma_R2

# =====================================================================
# KHỐI 2: GIỚI HẠN LỖI HELSTROM (CÔNG THỨC DẠNG ĐÓNG & LÝ TƯỞNG)
# =====================================================================
def helstrom_ideal(gamma_amp):
    """Xác suất lỗi Helstrom lý tưởng (không nhiễu) (Phương trình 41."""
    return 0.5 * (1 - np.sqrt(1 - np.exp(-4 * gamma_amp**2))) 

def helstrom_closed_form(gamma_amp, alpha, beta, sigma):
    """Xác suất lỗi Helstrom dạng đóng dưới fading Gamma-Gamma và nhiễu cộng (Phương trình 39)."""
    term_gamma2 = (4 * gamma_amp**2) / (1 + 8 * sigma**2) 
    ab = alpha * beta 
    
    coeff = 0.25 * np.sqrt(1 / (1 + 8 * sigma**2)) * (ab**beta) / \
            (gamma_func(alpha) * gamma_func(beta)) * (term_gamma2**(-alpha))
            
    part1 = (term_gamma2**(alpha - beta)) * gamma_func(alpha - beta) * gamma_func(beta) * \
            hyp1f1(beta, 1 - (alpha - beta), ab / term_gamma2) 
            
    part2 = (ab**(alpha - beta)) * gamma_func(beta - alpha) * gamma_func(alpha) * \
            hyp1f1(alpha, 1 + (alpha - beta), ab / term_gamma2) 
            
    return coeff * (part1 + part2)

def end_to_end_error(pe1, pe2):
    """Tổng xác suất lỗi đầu cuối của hệ thống tiếp sức 2 chặng (Phương trình 12)."""
    return pe1 * (1 - pe2) + (1 - pe1) * pe2 

# =====================================================================
# KHỐI 3: TOÁN TỬ MẬT ĐỘ TRONG CƠ SỞ FOCK VÀ HÀM MEIJER-G
# =====================================================================
def gamma_gamma_pdf_params(alpha, beta):
    """Hằng số phân phối Gamma-Gamma (Phương trình 10)[cite: 1]."""
    return 2 * (alpha * beta)**((alpha + beta) / 2) / (gamma_func(alpha) * gamma_func(beta)) 

def calc_rho_element(m, n, gamma_amp, alpha, beta, sigma, state_sign):
    """Tính phần tử (m, n) của ma trận mật độ thông qua hàm Meijer-G (Phương trình 24)[cite: 1]."""
    psi = state_sign * gamma_amp
    term1 = 1 / (sigma * np.sqrt(2 * np.pi)) 
    term2 = gamma_gamma_pdf_params(alpha, beta) 
    sigma_factor = 1 + 2 * sigma**2
    sum_val = mpmath.mpf(0)
    
    for k in range(m + n + 1):
        comb_coeff = comb(m + n, k) / np.sqrt(float(mpmath.fac(m) * mpmath.fac(n)))
        pow1 = (1 / sigma_factor)**((m + n - k - alpha - beta) / 2)
        pow2 = ((2 * sigma**2) / sigma_factor)**((k + 1) / 2)
        gamma_k = float(mpmath.gamma((k + 1) / 2)) 
        psi_pow = (abs(psi)**(-m - n - alpha - beta + k)) * (psi**(m + n - k)) 
        
        z = (alpha * beta * sigma_factor) / (psi**2) 
        a1 = 1 - (m + n - k + alpha + beta) / 2 
        b1, b2 = (alpha - beta) / 2, -(alpha - beta) / 2 
        
        meijer_g = mpmath.meijerg([[a1], []], [[b1, b2], []], z) 
        sum_val += comb_coeff * pow1 * pow2 * gamma_k * psi_pow * meijer_g 
        
    return float(term1 * term2 * sum_val) 

def helstrom_error_density_matrix(N, gamma_amp, alpha, beta, sigma):
    """Tính lỗi Helstrom bằng phân rã trị riêng ma trận mật độ cắt cụt (Phương trình 25, 26)[cite: 1]."""
    rho_0 = np.array([[calc_rho_element(m, n, gamma_amp, alpha, beta, sigma, 1) for n in range(N)] for m in range(N)]) 
    rho_1 = np.array([[calc_rho_element(m, n, gamma_amp, alpha, beta, sigma, -1) for n in range(N)] for m in range(N)]) 
    
    Delta = 0.5 * rho_1 - 0.5 * rho_0 
    eigenvalues = np.linalg.eigvalsh(Delta) 
    return 0.5 * (1 - np.sum(np.abs(eigenvalues))) 

# =====================================================================
# KHỐI 4: THỰC NGHIỆM MONTE CARLO LỚP VẬT LÝ
# =====================================================================
def monte_carlo_quantum_relay(gamma_amp, alpha, beta, sigma, num_bits=1000000):
    """Mô phỏng ngẫu nhiên kênh lượng tử với nhiễu Gamma-Gamma và Gaussian."""
    g1 = np.random.gamma(shape=alpha, scale=1.0/alpha, size=num_bits)
    g2 = np.random.gamma(shape=beta, scale=1.0/beta, size=num_bits)
    I_1 = g1 * g2
    
    epsilon_1 = np.random.normal(loc=0.0, scale=sigma, size=num_bits) 
    effective_amplitude = np.sqrt(I_1) * gamma_amp + epsilon_1 
    
    overlap = np.exp(-2 * (effective_amplitude**2))
    pe_array = 0.5 * (1 - np.sqrt(np.maximum(0, 1 - overlap**2))) 
    return np.mean(pe_array) 

# =====================================================================
# KHỐI 5: KỊCH BẢN THỰC THI KIỂM ĐỊNH TỔNG THỂ 
# =====================================================================
def run_all_simulations():
    print("Bắt đầu tái tạo kiểm định toàn bộ nghiên cứu\n")

    print("1/4. Đang tái tạo Kiến trúc Direct vs Relay (Hình 4)")
    sigma_noise_fig4 = 0.1
    gamma_vals_fig4 = np.linspace(1, 10, 20)
    direct_links = {2400: (3.3837, 1.5380), 1400: (3.9775, 1.8080)} 
    relay_links = {1200: (4.1658, 1.8935), 700: (4.8969, 2.2259)} 
    
    plt.figure(1, figsize=(10, 6))
    for i, L_direct in enumerate(direct_links.keys()):
        a_dir, b_dir = direct_links[L_direct]
        a_rel, b_rel = relay_links[L_direct // 2]
        
        pe_direct = [helstrom_closed_form(g, a_dir, b_dir, sigma_noise_fig4) for g in gamma_vals_fig4] 
        pe_relay_e2e = [end_to_end_error(helstrom_closed_form(g, a_rel, b_rel, sigma_noise_fig4), 
                                         helstrom_closed_form(g, a_rel, b_rel, sigma_noise_fig4)) for g in gamma_vals_fig4] 
        
        plt.semilogy(gamma_vals_fig4, pe_direct, marker='*', linestyle=':', label=f'Direct, L={L_direct}m')
        plt.semilogy(gamma_vals_fig4, pe_relay_e2e, marker='o', mfc='none', linestyle='-', label=f'Relay, L={L_direct//2}m')
    plt.title('Hình 4: Direct vs Quantum Relay')
    plt.grid(True, which="both", ls="--")
    plt.legend()

    # 2. Tái tạo Hình 6: Phân tích Bảo mật / Guessing Probability
    print("2/4. Đang tái tạo Phân tích rò rỉ bảo mật (Hình 6)")
    alpha_strong, beta_strong, sigma_noise_fig6 = 4.2, 1.4, 0.1 
    gamma_vals_fig6 = np.linspace(0.1, 3.0, 40)
    p_guess_ideal = [1 - helstrom_ideal(g) for g in gamma_vals_fig6] 
    p_guess_real = [1 - helstrom_closed_form(g, alpha_strong, beta_strong, sigma_noise_fig6) for g in gamma_vals_fig6] 
    
    plt.figure(2, figsize=(10, 6))
    plt.plot(gamma_vals_fig6, p_guess_ideal, 'b-^', label='Ideal (No Noise)')
    plt.plot(gamma_vals_fig6, p_guess_real, 'r-v', label='Real (Strong Turbulence)')
    plt.axhline(0.80, color='k', linestyle='--', label='80% Guessing Limit')
    plt.title('Hình 6: Relay Detection & Eavesdropper Confidence')
    plt.grid(True)
    plt.legend()

    # 3. Tái tạo Hình 7: Monte Carlo vs Analytical
    print("3/4. Đang chạy Monte Carlo lớp vật lý (Hình 7 - Weak Turbulence)")
    alpha_wk, beta_wk, sigma_noise_fig7 = 4.4, 2.0, 0.5 
    gamma_vals_fig7 = np.linspace(1.0, 5.0, 10)
    pe_ana_fig7 = [helstrom_closed_form(g, alpha_wk, beta_wk, sigma_noise_fig7) for g in gamma_vals_fig7]
    pe_sim_fig7 = [monte_carlo_quantum_relay(g, alpha_wk, beta_wk, sigma_noise_fig7, num_bits=500000) for g in gamma_vals_fig7]
    
    plt.figure(3, figsize=(10, 6))
    plt.semilogy(gamma_vals_fig7, pe_ana_fig7, 'b-', label='Analytical')
    plt.semilogy(gamma_vals_fig7, pe_sim_fig7, 'rs', mfc='none', label='Simulated')
    plt.title('Hình 7: Monte Carlo Verification (Weak Turbulence)')
    plt.grid(True, which="both", ls="--")
    plt.legend()

    # 4. Tái tạo tốc độ hội tụ Ma trận mật độ (Mô phỏng 1 điểm để tối ưu thời gian)
    print("4/4. Đang kiểm định tốc độ hội tụ không gian Fock (Hình 3)")
    alpha_mod, beta_mod, sigma_noise_fig3 = 4.0, 1.9, 1.0 
    gamma_test = 2.0
    pe_closed = helstrom_closed_form(gamma_test, alpha_mod, beta_mod, sigma_noise_fig3)
    print(f"Lỗi Helstrom (Dạng đóng Analytical): {pe_closed:.6f}")
    
    for N in [4, 8, 12]:
        pe_matrix = helstrom_error_density_matrix(N, gamma_test, alpha_mod, beta_mod, sigma_noise_fig3) 
        print(f"Lỗi Helstrom (Ma trận mật độ cắt cụt N={N}): {pe_matrix:.6f}")
        
    print("\nHoàn tất toàn bộ quy trình kiểm định.")
    plt.show()

if __name__ == "__main__":
    run_all_simulations()