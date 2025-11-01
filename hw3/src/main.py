# main.py
from runge_kutta import RungeKutta
from plotter import Plotter
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# 全局参数（用于 solve 和 solve_scipy）
P0 = 1e7
omega = np.pi / 100
t_max = 500
heatcoeff = 0.05


def cooling(t, T):
    T_val = T if T.ndim == 0 else T[0]
    if t < t_max:
        P = P0 * (1 + np.sin(omega * t))
    else:
        P = 0
    dTdt = -heatcoeff * (T_val - 25) + 0.001 * P
    return np.array([dTdt])


def rhs_scipy(t, T):
    if t < t_max:
        P = P0 * (1.0 + np.sin(omega * t))
    else:
        P = 0.0
    return -heatcoeff * (T - 25.0) + 0.001 * P


def solve():
    rk4 = RungeKutta(
        A=[[0, 0, 0, 0], [0.5, 0, 0, 0], [0, 0.5, 0, 0], [0, 0, 1, 0]],
        b=[1 / 6, 1 / 3, 1 / 3, 1 / 6],
        c=[0, 0.5, 0.5, 1],
        name="RK4",
    )
    t_vals, T_vals = rk4.solve(f=cooling, t_span=(0.0, 1000.0), y0=[50.0], h=1.0)
    T = T_vals[:, 0]
    print(f"[RK4] Max T: {T.max():.4f} °C at t={t_vals[T.argmax()]:.4f}s")

    for i in range(100, 1001, 100):
        print(f"t={t_vals[i]:.1f}s, T={T[i]:.4f}°C")

    return t_vals, T


def solve_scipy():
    sol = solve_ivp(
        rhs_scipy,
        (0, 1000),
        [50.0],
        method="RK45",
        t_eval=np.linspace(0, 1000, 1001),
        rtol=1e-8,
        atol=1e-10,
    )
    T = sol.y[0]
    print(f"[SciPy] Max T: {T.max():.4f} °C at t={sol.t[T.argmax()]:.4f}s")
    return sol.t, T


def plot_error_near_discontinuity(t_scipy, error, t_max=500):
    plt.figure(figsize=(8, 4))
    plt.semilogy(t_scipy, error, "g-", linewidth=1.2)
    plt.axvline(x=t_max, color="r", linestyle="--", label=r"$t_{\max} = 500$")
    plt.xlim(t_max - 50, t_max + 50)  # 490 to 510
    plt.xlabel("Time (s)")
    plt.ylabel("Absolute Error (°C)")
    plt.title("Error Spike Near $t = 500$ Due to Discontinuity in $P(t)$")
    plt.legend()
    plt.grid(True, linestyle=":", alpha=0.7)

    # 保存（确保目录存在）
    from pathlib import Path

    output_dir = Path("hw3/picture")
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / "error_spike.png", dpi=300)
    plt.close()


if __name__ == "__main__":
    # Solve both
    t_rk4, T_rk4 = solve()
    t_sp, T_sp = solve_scipy()

    # Plot comparison
    plotter = Plotter()
    plotter.plot_comparison(t_rk4, T_rk4, t_sp, T_sp, t_max)

    # Compute and plot error near discontinuity
    T_rk4_interp = np.interp(t_sp, t_rk4, T_rk4)
    error = np.abs(T_sp - T_rk4_interp)
    plot_error_near_discontinuity(t_sp, error, t_max)

    # Plot parameter sensitivity
    param_cases = [
        {"heatcoeff": 0.05, "P0": 1e7, "label": "P0=1e7, α=0.05"},
        {"heatcoeff": 0.02, "P0": 1e7, "label": "P0=1e7, α=0.02"},
        {"heatcoeff": 0.05, "P0": 2e7, "label": "P0=2e7, α=0.05"},
        {"heatcoeff": 0.02, "P0": 2e7, "label": "P0=2e7, α=0.02"},
    ]
    plotter.plot_param_sensitivity(
        param_sets=param_cases,
        t_span=(0, 1000),
        y0=[50.0],
        t_max=t_max,
        omega=omega,
        output_name="param_sensitivity.png",
    )
