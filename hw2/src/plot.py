import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

current_dir = os.path.dirname(os.path.abspath(__file__))
picture_dir = os.path.join(current_dir, "..", "picture")
os.makedirs(picture_dir, exist_ok=True)

x_data = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
y_data = np.array([3.5, 4.2, 3.8, 4.5, 4.0, 3.6])


cs_scipy = CubicSpline(x_data, y_data, bc_type="natural")

# ----------------------------
# a_i = y_i （i=0..5）
# b_i, c_i, d_i 对应区间 [x_i, x_{i+1}] 的系数
# S_i(x) = a_i + b_i*(x - x_i) + c_i*(x - x_i)^2 + d_i*(x - x_i)^3
# ----------------------------
a_coeffs = np.array([3.5, 4.2, 3.8, 4.5, 4.0, 3.6])
b_coeffs = np.array([0.054856459, -0.004712919, 0.008995215, 0.013732057, -0.033923445])
c_coeffs = np.array([0.0, -0.002978469, 0.003663876, -0.003427033, 0.001044258, 0.0])
d_coeffs = np.array(
    [-0.000049641, 0.000110706, -0.000118182, 0.000074522, -0.000017404]
)


def spline_manual(x_val):
    """
    使用你C++计算的系数进行插值
    """
    x_val = np.asarray(x_val)
    # 处理标量输入
    is_scalar = np.isscalar(x_val)
    if is_scalar:
        x_val = np.array([x_val])

    result = np.empty_like(x_val, dtype=float)

    for idx, x in enumerate(x_val):
        # 找到 x 所在的区间 i，使得 x_i <= x <= x_{i+1}
        if x < x_data[0] or x > x_data[-1]:
            raise ValueError(f"x={x} 超出插值区间 [{x_data[0]}, {x_data[-1]}]")
        i = np.searchsorted(x_data, x) - 1
        i = np.clip(i, 0, len(x_data) - 2)  # 确保 i 在 [0, 4]

        dx = x - x_data[i]
        result[idx] = (
            a_coeffs[i] + b_coeffs[i] * dx + c_coeffs[i] * dx**2 + d_coeffs[i] * dx**3
        )

    return result[0] if is_scalar else result


# ----------------------------
# 4. 生成密集点用于绘图和误差分析
# ----------------------------
x_fine = np.linspace(0, 100, 500)
y_scipy = cs_scipy(x_fine)
y_manual = spline_manual(x_fine)

# 计算误差
error = np.abs(y_scipy - y_manual)
max_error = np.max(error)

print(f"Maximum absolute error (scipy vs manual implementation): {max_error:.2e}")

# ----------------------------
# 5. 绘图
# ----------------------------
plt.figure(figsize=(10, 6))
plt.plot(
    x_data, y_data, "o", label="Original measurement point", markersize=8, color="red"
)
plt.plot(x_fine, y_scipy, "-", label="scipy", linewidth=2)
plt.plot(x_fine, y_manual, "--", label="C++", linewidth=2)

plt.xlabel("x (m)", fontsize=12)
plt.ylabel("Water Depth (m)", fontsize=12)
plt.title(
    "Comparison of Cubic Spline Interpolation: scipy vs manual (Natural Boundary)",
    fontsize=14,
)
plt.legend()
plt.ylim(3.0, 5.5)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig(
    os.path.join(picture_dir, "spline_comparison.png"), dpi=600, bbox_inches="tight"
)
plt.show()


plt.figure(figsize=(10, 4))
plt.plot(x_fine, error, "m-", linewidth=1.5)
plt.xlabel("x (m)", fontsize=12)
plt.ylabel("Absolute error", fontsize=12)
plt.title(
    f"Interpolation Error Distribution(≈ Maximum error {max_error:.2e})", fontsize=14
)
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig(os.path.join(picture_dir, "error.png"), dpi=600, bbox_inches="tight")
plt.show()
