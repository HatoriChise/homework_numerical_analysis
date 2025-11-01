import matplotlib.pyplot as plt
import numpy as np


class Ploter:
    @staticmethod
    def plot_convergence(history, title="Romberg Integration Convergence"):
        """
        绘制收敛曲线
        :param history: list of (iteration, integral_value)
        :param title: 图表标题
        """
        iters = [item[0] for item in history]
        values = [item[1] for item in history]

        plt.figure(figsize=(8, 5))
        plt.plot(iters, values, "bo-", label="Romberg Estimate")
        plt.xlabel("Iteration (k)")
        plt.ylabel("Integral Estimate")
        plt.title(title)
        plt.grid(True, linestyle="--", alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig("picture/convergence.png", dpi=300)
        plt.show()

    @staticmethod
    def plot_integrand():
        """
        f(x) = (4π²/5) * (1 - cos(πx/50)) * |sin(πx/50)| * sqrt(25 + 256*cos²(πx/50))
        """
        # 定义高分辨率 x 网格
        x = np.linspace(-100, 100, 4000)

        # 计算被积函数值（使用 abs 避免 sqrt(sin^2) 的数值问题）
        t = np.pi * x / 50
        f_x = (
            (4 * np.pi**2 / 5)
            * (1 - np.cos(t))
            * np.abs(np.sin(t))
            * np.sqrt(25 + 256 * np.cos(t) ** 2)
        )

        # 绘图
        plt.figure(figsize=(10, 5))
        plt.plot(x, f_x, "b-", linewidth=1.5)
        plt.xlabel("x", fontsize=12)
        plt.ylabel("f(x)", fontsize=12)
        plt.title("Integrand Function $f(x)$ on $[0, 100]$", fontsize=14)
        plt.grid(True, linestyle="--", alpha=0.6)
        plt.xlim(-100, 100)
        plt.ylim(bottom=0)
        plt.tight_layout()
        plt.savefig("hw1/picture/integrand.png", dpi=300)
        plt.show()
