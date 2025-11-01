import matplotlib.pyplot as plt


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
        plt.show()
