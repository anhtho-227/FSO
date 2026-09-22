"""
constants.py -- LỚP 0: Hằng số & Table I (Section IX.C, xác nhận bằng ảnh gốc trang 10)
"""

TABLE_I = {
    "weak":     (4.4, 2.0),
    "moderate": (4.0, 1.9),
    "strong":   (4.2, 1.4),
}

# Bảng α,β theo độ dài liên kết (Section IX.B, trang 9) -- đã kiểm tra tự nhất quán
DIRECT_LINKS = {2400: (3.3837, 1.5380), 1400: (3.9775, 1.8080),
                800: (4.7046, 2.1385), 400: (5.7921, 2.6328)}

RELAY_LINKS = {1200: (4.1658, 1.8935), 700: (4.8969, 2.2259),
               400: (5.7921, 2.6328), 200: (7.1309, 3.2413)}
