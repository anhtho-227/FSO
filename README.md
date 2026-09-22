# Tái tạo bài báo: Quantum Relay-Assisted Free-Space Optical Communication

Tái tạo đầy đủ mô hình toán học và các Hình 1, 3–7 của bài báo:

> Bhatnagar, Arti & Bhatnagar, "Quantum Relay-Assisted Free-Space Optical Communication",
> IEEE Photonics Journal, Vol. 17, No. 3, June 2025.

Trong quá trình tái tạo, đã phát hiện **4 lỗi toán học thật** trong bài báo — xác nhận
bằng ảnh gốc PDF (không dựa vào OCR text) kết hợp tính toán độc lập nhiều lớp. Chi tiết
đầy đủ (kèm ảnh chụp công thức gốc) xem tại [`analysis/ERRORS.md`](analysis/ERRORS.md).

## Cấu trúc thư mục

```
.
├── model/                    # Lõi toán học, 8 module = 8 lớp công thức (Section II–IX.D)
│   ├── constants.py             LỚP 0  Table I, bảng độ dài liên kết
│   ├── quantum_states.py        LỚP 1  Section II   — trạng thái lượng tử, Helstrom thuần khiết
│   ├── channel.py                LỚP 2  Section III  — Gamma-Gamma turbulence, Rytov variance
│   ├── relay.py                  LỚP 3  Section IV   — lỗi đầu-cuối relay 2 chặng
│   ├── density_operator.py      LỚP 4,5 Section V, App.A-B — toán tử mật độ (cố định & trung bình I1)
│   ├── helstrom.py               LỚP 6  Section VI-VII — Helstrom qua trị riêng
│   ├── closed_form.py            LỚP 7  Section VIII, App.C-E — công thức đóng xấp xỉ
│   └── security.py               LỚP 8  Section IX.D — Pe/Pguess bảo mật relay
│
├── analysis/                 # LỚP 9 — kiểm chứng & tài liệu lỗi
│   ├── self_tests.py             Bộ tự kiểm chứng (6 test, so khớp as-written vs corrected)
│   └── ERRORS.md                 Báo cáo đầy đủ 4 lỗi đã xác nhận (bằng ảnh gốc PDF)
│
├── data_generate/scripts/    # Sinh dữ liệu kênh dùng chung cho thí nghiệm
│   ├── generate_table1.py        Đối chiếu Table I với mô hình Rytov
│   └── generate_channel_samples.py   Sinh mẫu Gamma-Gamma + nhiễu Gauss
│
├── experiments/bpsk/         # LỚP 10 — mỗi Hình một thí nghiệm độc lập
│   ├── run_fig1.py
│   ├── run_fig3.py            (--mode match | --mode rigorous, xem bên dưới)
│   ├── run_fig4.py
│   ├── run_fig5.py
│   ├── run_fig6.py
│   └── run_fig7.py
│
├── results_final/            # Ảnh .png xuất ra sau khi chạy experiments/
├── run_all.py                 Chạy toàn bộ pipeline (self-test + mọi Hình)
└── .gitignore
```

## Cài đặt

```bash
pip install numpy scipy mpmath matplotlib
```

## Chạy

```bash
# Toàn bộ pipeline (kiểm chứng + tái tạo mọi Hình) — vài phút
python3 run_all.py

# Hoặc từng phần riêng:
python3 -m analysis.self_tests              # chỉ chạy bộ kiểm chứng (nhanh)
python3 experiments/bpsk/run_fig3.py --mode match      # khớp đúng độ lớn bài báo gốc
python3 experiments/bpsk/run_fig3.py --mode rigorous   # Helstrom THẬT, lộ rõ lỗi #4
```

## Tóm tắt 4 lỗi đã xác nhận (chi tiết đầy đủ + ảnh chụp trong `analysis/ERRORS.md`)

| # | Lỗi | Vị trí | Module tương ứng |
|---|---|---|---|
| 1 | Thiếu dấu trị tuyệt đối quanh Σλk | eq.(26)/(30) | `model/helstrom.py` |
| 2 | Thiếu điều kiện "k chẵn" | eq.(22)-(24)/(60) | `model/density_operator.py` |
| 3 | Thiếu hệ số ½ | eq.(22)-(24) | `model/density_operator.py` |
| 4 | Quy ước "nhiễu chia sẻ, đảo dấu" khiến Section VIII đánh giá thấp Pe thật 3–8 lần | eq.(32)-(39) | `model/closed_form.py` |

Mỗi công thức có lỗi trong `model/` đều có 2 phiên bản song song: hậu tố `_as_written`
(đúng y văn bản, tái lập lại đúng những gì bài báo in ra) và `_corrected` (đã sửa).

## Lưu ý quan trọng khi diễn giải kết quả

- `run_fig3.py --mode match` dùng cách tính "điều kiện hoá từng mẫu, nhiễu độc lập,
  Fock-basis cắt cụt cho trạng thái thuần khiết" — cách này **khớp đúng độ lớn** với
  Hình 3 gốc, nhưng **không phải** giới hạn Helstrom nghiêm ngặt cho toán tử mật độ
  đã trộn hoàn toàn.
- `run_fig3.py --mode rigorous` mới là giới hạn Helstrom **đúng chuẩn lý thuyết**
  (dùng `model/density_operator.py` + `model/helstrom.py` đã sửa lỗi #2,#3,#1) —
  cho kết quả **cao hơn 3–8 lần** so với công thức đóng eq.(39) của bài báo.
- Hình 4, 5, 6, 7 hiện dùng eq.(39) (`Pe1_closed_form_eq39`) để nhất quán với cách
  bài báo trình bày — con số tuyệt đối trong các hình này do đó **thấp hơn giới hạn
  Helstrom thật ~3-8 lần** (kế thừa lỗi #4), dù xu hướng định tính nhìn chung vẫn đúng.
