import numpy as np
from integrator import Integrator
from plotter import Ploter


def integrand(x):
    """被积函数"""
    return (
        (4 * np.pi**2 / 5)
        * (1 - np.cos(np.pi * x / 50))
        * np.abs(np.sin(np.pi * x / 50))
        * np.sqrt(25 + 256 * np.cos(np.pi * x / 50) ** 2)
    )


def main():
    a, b = 0, 100
    tol = 1e-12
    max_iter = 10  # Romberg 表最大 8x8（通常 5~7 次已足够高精度）

    # 创建积分器
    integrator = Integrator(func=integrand, a=a, b=b, tol=tol, max_iter=max_iter)

    # 执行 Romberg 积分，获取收敛历史
    history = integrator.romberg()

    # 打印结果
    print("Romberg Integration History:")
    print(f"{'Iter':<6} {'Estimate':<20}")
    print("-" * 30)
    for it, val in history:
        print(f"{it:<6} {val:<20.12f}")

    final_value = history[-1][1]
    print(f"\nFinal Romberg Estimate: {final_value:.12f}")

    # 绘制收敛曲线
    Ploter.plot_convergence(
        history, title="Romberg Integration Convergence for Given Integral"
    )

    Ploter.plot_integrand()


if __name__ == "__main__":
    main()
