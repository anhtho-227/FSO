"""
run_all.py -- Chạy TOÀN BỘ pipeline tái tạo: bộ tự kiểm chứng + mọi thí nghiệm/hình.
Gọi từ thư mục gốc: python3 run_all.py
"""
import subprocess
import sys
import os

ROOT = os.path.dirname(os.path.abspath(__file__))


def run(cmd):
    print(f"\n>>> {' '.join(cmd)}")
    subprocess.run(cmd, cwd=ROOT, check=True)


def main():
    run([sys.executable, "-m", "analysis.self_tests"])
    run([sys.executable, "experiments/bpsk/run_fig1.py"])
    run([sys.executable, "experiments/bpsk/run_fig3.py", "--mode", "match"])
    run([sys.executable, "experiments/bpsk/run_fig3.py", "--mode", "rigorous"])
    run([sys.executable, "experiments/bpsk/run_fig4.py"])
    run([sys.executable, "experiments/bpsk/run_fig5.py"])
    run([sys.executable, "experiments/bpsk/run_fig6.py"])
    run([sys.executable, "experiments/bpsk/run_fig7.py"])
    print("\n=== HOAN TAT. Xem ket qua trong results_final/ ===")


if __name__ == "__main__":
    main()
