import numpy as np
from typing import Callable, Tuple, Optional


class RungeKutta:
    def __init__(self, A, b, c, name: str = "RK"):
        self.A = np.array(A, dtype=float)
        self.b = np.array(b, dtype=float)
        self.c = np.array(c, dtype=float)
        self.s = len(self.b)
        self.name = name

        assert self.A.ndim == 2 and self.A.shape == (self.s, self.s), "A must be s×s"
        assert len(self.c) == self.s, "Length of c must be s"
        assert np.allclose(self.c[0], 0.0), "c[0] must be 0"
        if not np.allclose(np.triu(self.A), 0):
            raise ValueError("Only explicit Runge-Kutta methods are supported.")

    def solve(
        self,
        f: Callable[[float, np.ndarray], np.ndarray],
        t_span: Tuple[float, float],
        y0,
        h: Optional[float] = None,
        n_steps: Optional[int] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        if h is None and n_steps is None:
            raise ValueError("Either 'h' or 'n_steps' must be specified.")
        if h is not None and n_steps is not None:
            raise ValueError("Specify only one of 'h' or 'n_steps'.")

        t0, tf = t_span
        y = np.atleast_1d(np.asarray(y0, dtype=float))
        t = float(t0)

        if h is None:
            h = (tf - t0) / n_steps

        t_vals = [t]
        y_vals = [y.copy()]

        while t < tf:
            current_h = min(h, tf - t)
            k = []
            for i in range(self.s):
                if i == 0:
                    ki = f(t, y)
                else:
                    dy = np.zeros_like(y)
                    for j in range(i):
                        dy += self.A[i, j] * k[j]
                    ti = t + self.c[i] * current_h
                    yi = y + current_h * dy
                    ki = f(ti, yi)
                k.append(ki)

            y += current_h * sum(self.b[i] * k[i] for i in range(self.s))
            t += current_h
            t_vals.append(t)
            y_vals.append(y.copy())

        return np.array(t_vals), np.array(y_vals)
