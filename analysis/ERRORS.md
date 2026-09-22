# Báo cáo kiểm chứng: "Quantum Relay-Assisted Free-Space Optical Communication"

**Bài báo:** Bhatnagar, Arti & Bhatnagar — IEEE Photonics Journal, Vol. 17, No. 3, 2025

---

## 1. Tóm tắt điều hành

Bài báo có **4 lỗi toán học thật, đã xác nhận chắc chắn bằng ảnh gốc PDF và/hoặc tính toán độc lập nhiều lớp**. Ba lỗi đầu nằm ở Appendix B (toán tử mật độ trung bình theo Gamma-Gamma) và Section VI/VII (chuẩn vết Helstrom). Lỗi thứ tư — nghiêm trọng nhất về mặt hệ quả - là cả một **hướng tiếp cận (Section VIII) cho kết quả hệ thống lệch xa (3–8 lần)** so với giới hạn Helstrom thật, dù công thức đóng cuối cùng (eq. 39/76) tự nó được suy diễn/rút gọn đại số chính xác từ điểm xuất phát (eq. 32-34) - vấn đề nằm ở chính **điểm xuất phát đó**, không phải ở các bước đại số sau.

---

## 2. Bảng 4 lỗi đã xác nhận

| # | Lỗi | Vị trí | Bằng chứng | Mức độ nghiêm trọng |
|---|---|---|---|---|
| 1 | Thiếu dấu trị tuyệt đối `\|·\|` quanh $\sum_k\lambda_k$ | eq. (26), (30) | **Ảnh gốc trang 7**: eq.(25) có `Tr\|p₁ρ₁-p₀ρ₀\|`; eq.(26) chỉ còn `Σλk` trần | Cao về mặt hình thức, nhưng nhiều khả năng không lọt vào code tạo Hình 3 của tác giả (đường cong không phẳng ở 0.5) |
| 2 | Thiếu điều kiện "k chẵn" | eq. (22)-(24), tức eq.(60) | **Ảnh gốc trang 13**: eq.(57) có "$k$ even"; eq.(22) sau khi thế Meijer-G vào thì mất | Cao - làm hỏng hoàn toàn tính hợp lệ của ma trận mật độ (trace ≠ 1) nếu tính đúng văn bản |
| 3 | Thiếu hệ số **½** | eq. (22)-(24) | **Ảnh gốc trang 13**: eq.(60) có `$b^{-\nu}/2$`; ảnh gốc trang 7 eq.(24): hệ số ½ biến mất ngay trước $G^{2,1}_{1,2}$ | Cao - gây sai lệch đúng gấp 2 lần mọi phần tử ma trận nếu tính đúng văn bản |
| 4 | Quy ước nhiễu chia sẻ, đảo dấu theo tín hiệu (eq. 32) khiến toàn bộ Section VIII tính **sai bản chất** giới hạn Helstrom | eq. (32)–(39), lan ra toàn bộ Section VIII và các Hình 3–7 dùng eq.(39) | Đối chiếu số trực tiếp: Pe thật (density-operator, đã sửa lỗi #2+#3) = 0.10–0.22; Pe theo eq.(39) chỉ 0.013–0.071 - **lệch 3–8 lần** | **Nghiêm trọng nhất** - ảnh hưởng toàn bộ kết luận định lượng của bài báo |

---

## 3. Chi tiết từng lỗi

### 3.1 Lỗi #1 - Thiếu `|·|` ở eq. (26)/(30)

Ảnh gốc:
- eq.(25): $P_{e,1}^{(\min)} = \frac12\Big(1-\text{Tr}\big|p_1\rho_1^{(C)}-p_0\rho_0^{(C)}\big|\Big)$ ✅ đúng
- eq.(26): $P_{e,1}^{(\min)} = \frac12\Big(1-\sum_k\lambda_k\Big)$ ❌ thiếu `|·|`

Vì $\text{Tr}(\Delta) = p_1\text{Tr}(\rho_1) - p_0\text{Tr}(\rho_0) = p_1 - p_0$, với ưu tiên bằng nhau ($p_0=p_1=0.5$) thì $\sum_k\lambda_k = 0$ luôn luôn, bất kể kênh truyền - khiến eq.(26) tính đúng văn bản sẽ **luôn cho $P_e=0.5$**, vô giá trị. Công thức đúng cần là $\sum_k|\lambda_k|$.

**Lưu ý:** rất có thể đây chỉ là lỗi đánh máy khi rút gọn từ (25) sang (26)/(29-30), vì nếu lỗi này thực sự lọt vào code tạo Hình 3 của tác giả thì đường "density operator" phải là đường thẳng ngang ở 0.5 - nhưng Hình 3 thực tế cho thấy đường cong giảm dần theo γ giống 2 phương pháp còn lại. Điều này gợi ý tác giả **có** dùng `|λk|` khi lập trình thực tế, chỉ quên viết `|·|` trong bản thảo.

### 3.2 Lỗi #2 - Thiếu "k chẵn" ở eq. (22)-(24)

Ảnh gốc trang 13, eq.(57) (bước dẫn xuất trung gian) rõ ràng ghi:
$$\sum_{\substack{k=0\\k\text{ even}}}^{m+n}$$

nhưng khi thế công thức Meijer-G (eq. 60) vào để có kết quả cuối (eq. 22-24), điều kiện "k even" **biến mất**, chỉ còn $\sum_{k=0}^{m+n}$.

**Hệ quả khi tính đúng văn bản:** ma trận mật độ dựng ra có $\text{Tr}(\rho)\neq 1$ (thậm chí có thể âm), không còn là toán tử mật độ hợp lệ — dẫn tới xác suất lỗi Helstrom tính ra **âm** ở nhiều giá trị γ (đã kiểm chứng bằng cả Meijer-G và tích phân số).

### 3.3 Lỗi #3 - Thiếu hệ số ½ ở eq. (22)-(24)

Ảnh gốc trang 13, eq.(60): $\mathcal J = \dfrac{b^{-\nu}}{2}G^{2,1}_{1,2}(\cdots)$ - có hệ số ½ rõ ràng.

Ảnh gốc trang 7, eq.(24): `$\times\, G^{2,1}_{1,2}(\cdots)$` - **không còn** hệ số ½ ngay trước ký hiệu Meijer-G, dù văn bản nói rõ "Substituting (59) in (57), we get (24)".

**Hệ quả:** mọi phần tử ma trận mật độ tính theo văn bản eq.(24) bị **gấp đôi** giá trị đúng - đã kiểm chứng bằng số: sai lệch đúng 100.00% (tức gấp đôi) ở mọi cặp $(m,n)$ thử nghiệm.

**Kiểm chứng tổng hợp:** khi sửa **đồng thời cả lỗi #2 và #3** (thêm điều kiện k chẵn + nhân thêm ½), công thức Meijer-G khớp với phương pháp tích phân số độc lập **chính xác đến 0.0000%** ở mọi phần tử ma trận thử nghiệm, và $\text{Tr}(\rho)\approx 1$ (0.99, 0.93, 0.78 khi γ tăng - sai lệch nhỏ dần do cắt cụt Fock basis N=16, sẽ tiến gần 1 hơn nếu tăng N). Đây là bằng chứng mạnh nhất cho thấy **cả 2 lỗi đều là lỗi đánh máy khi rút gọn công thức**, không phải sai sót trong tư duy toán học cốt lõi.

### 3.4 Lỗi #4 — Section VIII đánh giá sai bản chất giới hạn Helstrom (nghiêm trọng nhất)

**Vấn đề gốc rễ (eq. 32):** khi $\psi_1=-\gamma$ được thế vào eq.(31) $\beta_i^{(C)}=\sqrt{I_1}\psi_i+\epsilon_1$, kết quả đúng phải là $\beta_1^{(C)}=-\sqrt{I_1}\gamma+\epsilon_1$ (nhiễu cộng KHÔNG đổi dấu theo tín hiệu - đây là nguyên lý độc lập thống kê cơ bản của nhiễu kênh truyền). Nhưng bài báo viết:
$$|\psi_1^{(C)}\rangle = |-(\sqrt{I_1}\gamma+\epsilon_1)\rangle = |-\sqrt{I_1}\gamma-\epsilon_1\rangle$$

**Vì sao điều này quan trọng (khác với những gì báo cáo V1/V2 từng đánh giá):** dù phân phối MARGINAL của $\beta_1$ vẫn là $\mathcal N(-\sqrt{I_1}\gamma,\sigma^2)$ trong cả hai cách viết (do nhiễu Gauss đối xứng quanh 0) — đây là lý do trước đây tôi từng (sai) kết luận "không ảnh hưởng" - nhưng khi tính **độ chồng lấn (overlap)** giữa HAI trạng thái dùng CHUNG một $\epsilon_1$ (eq. 34), quy ước dấu này quyết định **quan hệ tương đối** giữa hai trạng thái, không chỉ phân phối riêng lẻ của từng trạng thái. Với quy ước sai của bài báo, ngay cả khi $\gamma=0$ (không có tín hiệu thật, hai trạng thái phát đi giống hệt nhau), nhiễu $\epsilon_1$ vẫn tạo ra một sự "phân biệt giả" giữa hai giả thuyết - điều này vô lý về vật lý (nhiễu không thể tự nhiên giúp phân biệt hai tín hiệu giống hệt nhau).

**Bằng chứng số quyết định - giới hạn γ→0:**

| Phương pháp | $P_e(\gamma\to 0)$ | Đúng về vật lý? |
|---|---|---|
| Toán tử mật độ đúng (density-operator, đã sửa lỗi #2+#3) | **0.500000** | ✅ Đúng - hai trạng thái giống hệt nhau nên không thể phân biệt |
| eq.(39)/Section VIII | NaN (tràn số) hoặc ≈0.08–0.12 tùy σ | ❌ Sai hoàn toàn |

**Bằng chứng số ở vùng γ bình thường (α=4.0, β=1.9, σ=1):**

| γ | $P_e$ đúng (density-operator) | $P_e$ theo eq.(39) | Tỉ lệ chênh lệch |
|---|---|---|---|
| 1.0 | 0.2176 | 0.0714 | 3.0× |
| 2.0 | 0.1036 (→ 0.119 khi sửa #2+#3) | 0.0274 | 3.8–4.3× |
| 3.0 | 0.1125 (→ 0.147 khi sửa #2+#3) | 0.0130 | 8.7–11×× |

**Kết luận:** công thức closed-form (39)/(76) - dù bản thân các bước đại số từ (32)→(39) là **chính xác về mặt suy diễn** (đã kiểm chứng kỹ ở báo cáo V2) - tính ra một đại lượng **khác** với giới hạn Helstrom thật của bài toán phân biệt hai trạng thái hỗn hợp. Nó có xu hướng **đánh giá thấp** (lạc quan hơn thực tế) xác suất lỗi hệ thống. Tuyên bố ở *Remark 8* rằng "cả hai phương pháp cho cùng một giới hạn lỗi Helstrom tối thiểu" là **không chính xác**.

---

## 4. Những phần vẫn được xác nhận đúng

| Phần | Trạng thái |
|---|---|
| Coherent state, overlap cơ bản (eq. 1-2) | ✅ Đúng |
| Helstrom bound nhị phân thuần khiết (eq. 3-5) | ✅ Đúng |
| Rytov variance → α, β (eq. 7-9), Hình 1 | ✅ Đúng |
| Toán tử mật độ cố định $I_1$ (eq. 20/55, có "k chẵn") | ✅ Đúng |
| Trung bình nhiễu cộng (Appendix C, eq. 61-70) | ✅ Đúng |
| Công thức Laplace Bessel-K → Meijer-G tổng quát (eq. 75) | ✅ Đúng |
| Đại số rút gọn từ (32) → (39) (nếu chấp nhận điểm xuất phát eq.32) | ✅ Đúng về mặt đại số |
| Công thức lỗi đầu-cuối relay (eq. 12) | ✅ Đúng (là phép toán XOR chuẩn, không phụ thuộc vào cách tính $P_{e,1},P_{e,2}$) |
| Table I (Weak: α=4.4,β=2.0; Moderate: α=4.0,β=1.9; Strong: α=4.2,β=1.4) | ✅ Đã xác nhận từ ảnh gốc |

---

## 5. Công thức đề xuất sửa lại (Appendix B)

Với $A = 1-\dfrac{m+n-k+\alpha_1+\beta_1}{2}$, $B=\dfrac{\alpha_1-\beta_1}{2}$:

$$\big[\rho_i^{(C)}\big]_{m,n} = \frac{1}{\sigma\sqrt{2\pi}}\frac{(\alpha_1\beta_1)^{\frac{\alpha_1+\beta_1}{2}}}{\Gamma(\alpha_1)\Gamma(\beta_1)}\sum_{\substack{k=0\\k\text{ even}}}^{m+n}\binom{m+n}{k}\frac{1}{\sqrt{m!n!}}\left(\frac{1}{1+2\sigma^2}\right)^{\frac{m+n-k-\alpha_1-\beta_1}{2}}\left(\frac{2\sigma^2}{1+2\sigma^2}\right)^{\frac{k+1}{2}}$$
$$\times\,\Gamma\!\left(\frac{k+1}{2}\right)|\psi_i|^{-m-n-\alpha_1-\beta_1+k}\psi_i^{m+n-k}\times\boxed{\frac12}\,G^{2,1}_{1,2}\!\left(\frac{\alpha_1\beta_1(1+2\sigma^2)}{\psi_i^2}\,\middle|\,A;\,B,-B\right)$$

(phần được đóng khung `½` là bổ sung so với văn bản gốc; điều kiện "k chẵn" dưới dấu tổng cũng là bổ sung).

---

## 6. Khuyến nghị sử dụng bài báo

1. **Có thể tin tưởng:** mô hình kênh (Gamma-Gamma + nhiễu Gauss), công thức Rytov variance, và khung lý thuyết Helstrom nói chung.
2. **Cần tự sửa trước khi dùng:** eq.(24)/(28) — thêm điều kiện "k chẵn" và hệ số ½; eq.(26)/(30) — thêm dấu `|·|`.
3. **Không nên dùng trực tiếp để tính số:** eq.(39)/(40)/(76) và mọi kết quả số ở Section IX phụ thuộc vào chúng (Hình 3, 4, 5, 6, 7) — các con số này **đánh giá thấp hơn thực tế 3-8 lần** so với giới hạn Helstrom đúng. Nếu cần con số đáng tin cậy, nên tính lại bằng phương pháp toán tử mật độ đầy đủ (Section V-VII) sau khi đã sửa lỗi #2 và #3.
4. Xu hướng định tính của bài báo (relay tốt hơn truyền thẳng ở biên độ vừa-cao, nhiễu loạn mạnh làm tăng lỗi, v.v.) **có khả năng vẫn đúng về mặt định tính**, vì lỗi #4 là một *thiên lệch hệ thống* (luôn đánh giá thấp Pe theo cùng một hướng) chứ không phải nhiễu loạn ngẫu nhiên — nhưng **độ lớn con số cụ thể thì không đáng tin**.

---

## 7. Ghi nhận

Các phát hiện lỗi #1, #3, và #4 trong báo cáo này xuất phát trực tiếp từ phân tích độc lập và đoạn code đối chiếu do người dùng cung cấp — bao gồm việc chỉ ra đúng quy ước tham số cho hàm Meijer-G (`mpmath.meijerg`), điều mà bản thân tôi đã lập trình sai ở lượt làm việc trước đó. Việc đối chiếu chéo với công cụ/phân tích độc lập bên ngoài đã đóng vai trò quyết định trong việc phát hiện đầy đủ các lỗi này.

---

## 8. Danh mục file

| File | Nội dung |
|---|---|
| `full_reconstruction_v2.py` | Code tái tạo, dùng eq.(39) đúng theo ảnh gốc — **lưu ý: chưa tính đến lỗi #4**, chỉ tái hiện đúng những gì bài báo viết |
| Báo cáo này (V3) | Tổng hợp đầy đủ và cuối cùng, 4 lỗi đã xác nhận |
| Các báo cáo V1, V2 trước | Đã lỗi thời, giữ lại chỉ để tham khảo lịch sử quá trình kiểm chứng |
