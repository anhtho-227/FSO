"""
model -- Bộ mô hình toán học tái tạo từ bài báo "Quantum Relay-Assisted FSO
Communication" (Bhatnagar et al., IEEE Photonics Journal, 2025).

10 lớp, mỗi module tương ứng 1-2 lớp theo đúng cấu trúc Section của bài báo:
    quantum_states.py    LỚP 1  (Section II)
    channel.py           LỚP 2  (Section III)
    relay.py             LỚP 3  (Section IV)
    density_operator.py  LỚP 4, 5  (Section V, Appendix A-B)
    helstrom.py          LỚP 6  (Section VI-VII)
    closed_form.py       LỚP 7  (Section VIII, Appendix C-E)
    security.py          LỚP 8  (Section IX.D)

Xem ERRORS.md trong thư mục analysis/ để biết đầy đủ 4 lỗi đã xác nhận.
"""
from .constants import TABLE_I, DIRECT_LINKS, RELAY_LINKS
from .quantum_states import coherent_overlap, helstrom_binary_pure
from .channel import rytov_variance, alpha_beta_from_rytov, gamma_gamma_pdf, sample_gamma_gamma
from .relay import end_to_end_error
from .density_operator import (
    rho_element_fixed_I,
    rho_element_numeric_reference,
    rho_element_meijerG_as_written,
    rho_element_meijerG_corrected,
    build_rho_matrix,
)
from .helstrom import helstrom_bound_as_written, helstrom_bound_corrected
from .closed_form import (
    overlap_shared_noise_as_written,
    Pe_exact_conditional_eq36,
    Pe_approx_conditional_eq37,
    Pe1_averaged_over_noise_eq38,
    Pe1_closed_form_eq39,
    Pe1_closed_form_hp,
    Pe1_closed_form_robust,
)
from .security import Pe_ideal, Pguess_ideal, Pe_real, Pguess_real
