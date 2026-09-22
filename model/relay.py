"""
relay.py -- LỚP 3: Hệ thống relay 2 chặng (Section IV)  [ĐÃ XÁC NHẬN ĐÚNG]
"""


def end_to_end_error(pe1, pe2):
    """P_e = Pe1(1-Pe2) + (1-Pe1)Pe2   -- eq.(12)
    Công thức XOR chuẩn -- đúng bất kể Pe1, Pe2 được tính bằng phương pháp nào."""
    return pe1*(1-pe2) + (1-pe1)*pe2
