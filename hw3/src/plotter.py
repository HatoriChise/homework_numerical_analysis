# plotter.py
import matplotlib

matplotlib.use("Agg")  # Use non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp


class Plotter:
    def __init__(self, output_dir: str = "hw3/picture"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def plot_comparison(self, t_rk4, T_rk4, t_scipy, T_scipy, t_max: float):
        """
        Plot comparison between custom RK4 and SciPy RK45 solutions.
        """
        fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Temperature comparison
        axs[0].plot(t_rk4, T_rk4, "b-", linewidth=1.2, label="Custom RK4 (fixed h=1.0)")
        axs[0].plot(
            t_scipy, T_scipy, "r--", linewidth=1.2, label="SciPy RK45 (adaptive)"
        )
        axs[0].axvline(
            x=t_max, color="k", linestyle=":", linewidth=1.0, label=r"$t_{\max}$"
        )
        axs[0].set_ylabel("Temperature (°C)")
        axs[0].set_title("Comparison of RK4 and SciPy RK45 Solutions")
        axs[0].legend()
        axs[0].grid(True, linestyle=":", alpha=0.7)

        # Absolute error
        # Interpolate RK4 onto SciPy time points for fair comparison
        T_rk4_interp = np.interp(t_scipy, t_rk4, T_rk4)
        error = np.abs(T_scipy - T_rk4_interp)
        axs[1].plot(t_scipy, error, "g-", linewidth=1.0)
        axs[1].axvline(x=t_max, color="k", linestyle=":", linewidth=1.0)
        axs[1].set_xlabel("Time (s)")
        axs[1].set_ylabel("Absolute Error (°C)")
        axs[1].set_yscale("log")
        axs[1].set_title("Absolute Error (SciPy - RK4)")
        axs[1].grid(True, linestyle=":", alpha=0.7)

        plt.tight_layout()
        plt.savefig(self.output_dir / "comparison.png", dpi=300)
        plt.close()

    def plot_param_sensitivity(
        self, param_sets, t_span, y0, t_max, omega, output_name="sensitivity.png"
    ):
        """
        Plot temperature response under different (heatcoeff, P0) settings.

        param_sets: list of dict, e.g.,
            [{"heatcoeff": 0.05, "P0": 1e4, "label": "Case 1"}, ...]
        """
        fig, ax = plt.subplots(figsize=(10, 6))

        for params in param_sets:
            heatcoeff = params["heatcoeff"]
            P0 = params["P0"]
            label = params.get("label", f"P0={P0:.0e}, α={heatcoeff}")

            def rhs(t, T):
                if t < t_max:
                    P = P0 * (1.0 + np.sin(omega * t))
                else:
                    P = 0.0
                return -heatcoeff * (T - 25.0) + 0.001 * P

            sol = solve_ivp(
                rhs,
                t_span,
                y0,
                method="RK45",
                t_eval=np.linspace(t_span[0], t_span[1], 1001),
                rtol=1e-8,
                atol=1e-10,
            )
            ax.plot(sol.t, sol.y[0], label=label)

        ax.axvline(
            x=t_max, color="k", linestyle=":", linewidth=1.0, label=r"$t_{\max}$"
        )
        ax.set_xlabel("Time (s)")
        ax.set_ylabel("Temperature (°C)")
        ax.set_title("Sensitivity to Heat Transfer Coefficient and Power Amplitude")
        ax.legend()
        ax.grid(True, linestyle=":", alpha=0.7)
        plt.tight_layout()
        plt.savefig(self.output_dir / output_name, dpi=300)
        plt.close()
