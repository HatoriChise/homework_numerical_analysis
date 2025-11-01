import numpy as np


class Integrator:
    def __init__(self, func, a, b, tol=1e-10, max_iter=10):
        """
        initalize the integrator
        func: integrand function
        a: lower bound
        b: upper bound
        tol: error tolerance
        max_iter: maximum number of iterations
        """
        self.func = func
        self.a = a
        self.b = b
        self.tol = tol
        self.max_iter = max_iter

    def romberg(self):
        """
        compute the integral using Romberg integration
        return: history of iterations and results
        """
        R = [[0.0 for _ in range(self.max_iter)] for _ in range(self.max_iter)]
        history = []

        # 第0列：复合梯形公式
        h = self.b - self.a
        R[0][0] = 0.5 * h * (self.func(self.a) + self.func(self.b))
        history.append((1, R[0][0]))

        for k in range(1, self.max_iter):
            # 计算 T_k（第k行第0列）
            h /= 2.0
            sum_f = 0.0
            for i in range(1, 2**k, 2):
                x = self.a + i * h
                sum_f += self.func(x)
            R[k][0] = 0.5 * R[k - 1][0] + h * sum_f

            # Richardson 外推：填充第k行
            for j in range(1, k + 1):
                factor = 4**j
                R[k][j] = (factor * R[k][j - 1] - R[k - 1][j - 1]) / (factor - 1)

            history.append((k + 1, R[k][k]))

            # 检查收敛：比较 R[k][k] 与 R[k-1][k-1]
            if k >= 1:
                if abs(R[k][k] - R[k - 1][k - 1]) < self.tol:
                    break

        return history
