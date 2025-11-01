import numpy as np
from scipy.integrate import quad


# 定义被积函数 f(x)
def integrand(x):
    term1 = 20 * np.pi
    term2 = 1 - np.cos(np.pi * x / 50)

    sq_term_a = (np.pi / 5 * np.sin(np.pi * x / 50)) ** 2
    sq_term_b = (-8 * np.pi / 25 * np.sin(np.pi * x / 25)) ** 2

    sqrt_term = np.sqrt(sq_term_a + sq_term_b)

    return term1 * term2 * sqrt_term


# 积分的下限和上限
a = 0
b = 100

# 使用scipy.integrate.quad进行数值积分
# quad返回一个包含 (积分值, 估计误差) 的元组
integral_value, error_estimate = quad(integrand, a, b)

# 打印结果
print(f"积分值: {integral_value}")
print(f"估计误差: {error_estimate}")

# 例如，保留三位小数的格式化输出
print(f"积分值的近似值 (3位小数): {integral_value:.6f}")
