import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu_factor, lu_solve


class SystemEquation:
    def __init__(self, A, b):
        """
        Initialize the system Ax = b.

        Parameters:
        A (np.ndarray): Coefficient matrix (n x n)
        b (np.ndarray): Right-hand side vector (n,)
        """
        self.A = A.copy()
        self.b = b.copy()
        self.n = A.shape[0]

        if A.shape[0] != A.shape[1]:
            raise ValueError("Matrix A must be square.")
        if b.shape[0] != self.n:
            raise ValueError("Vector b length must match matrix size.")

    def lu_solve(self):
        """
        Solve Ax = b using LU decomposition (direct method).

        Returns:
        np.ndarray: Solution vector x
        """
        lu, piv = lu_factor(self.A)
        x = lu_solve((lu, piv), self.b)
        return x

    def sor_solve(self, omega=1.0, tol=1e-8, max_iter=1000):
        """
        Solve Ax = b using Successive Over-Relaxation (SOR) method.

        Parameters:
        omega (float): Relaxation parameter (0 < omega < 2)
        tol (float): Tolerance for convergence (residual norm)
        max_iter (int): Maximum number of iterations

        Returns:
        tuple: (x, residuals, solutions)
               x: Final solution
               residuals: List of residual norms at each iteration
               solutions: List of solution vectors at each iteration
        """
        if not 0 < omega < 2:
            raise ValueError("Relaxation parameter omega must be in (0, 2).")

        x = np.zeros(self.n)  # Initial guess
        residuals = [np.linalg.norm(self.A @ x - self.b)]  # Initial residual
        solutions = [x.copy()]  # Initial solution

        for iter in range(max_iter):
            x_old = x.copy()

            # SOR iteration
            for i in range(self.n):
                sigma = np.dot(self.A[i, :i], x[:i]) + np.dot(
                    self.A[i, i + 1 :], x_old[i + 1 :]
                )
                x[i] = (1 - omega) * x_old[i] + (omega / self.A[i, i]) * (
                    self.b[i] - sigma
                )

            # Store residual and solution AFTER update
            residual = np.linalg.norm(self.A @ x - self.b)
            residuals.append(residual)
            solutions.append(x.copy())

            # Check convergence (use strict inequality)
            if residual < tol:
                break

        return x, residuals, solutions

    def plot_convergence(
        self, residuals, title="Convergence Plot", exact_solution=None
    ):
        """
        Plot the convergence history (residual norm vs iterations).

        Parameters:
        residuals (list): List of residual norms
        title (str): Plot title
        exact_solution (np.ndarray, optional): Exact solution for error plot
        """
        plt.figure(figsize=(8, 5))
        plt.semilogy(
            range(len(residuals)), residuals, "b-", linewidth=2, label="Residual Norm"
        )

        if exact_solution is not None:
            errors = [np.linalg.norm(x - exact_solution) for x in solutions]
            plt.semilogy(
                range(len(errors)), errors, "r--", linewidth=2, label="Error Norm"
            )

        plt.xlabel("Iteration")
        plt.ylabel("Norm (log scale)")
        plt.title(title)
        plt.grid(True, which="both", linestyle="--")
        plt.legend()
        plt.show()


def main():
    # Define A and b (from your problem)
    A = np.array(
        [
            [1.05, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0.007, -0.015, 0, 0, 0, 0, 0, 0, 0, 0],
            [-0.9, 0, 1, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0.017, -0.04, 0, 0, 0, 0, 0, 0],
            [0, 0, -0.85, 0, 1, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 0.04, -0.103, 0, 0, 0, 0],
            [0, 0, 0, 0, -0.92, 0, 1, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0.105, -0.25, 0, 0],
            [0, 0, 0, 0, 0, 0, -0.88, 0, 1, 1],
            [0, 0, 0, 0, 0, 0, 0, 0, 0.25, -0.808],
        ],
        dtype=float,
    )
    b = np.zeros(10)
    b[0] = 500.0

    # Create a SystemEquation object
    system = SystemEquation(A, b)

    # Solve using LU decomposition
    x_lu = system.lu_solve()
    print("LU Solution:", x_lu)

    # Solve using SOR method
    x_sor, residuals, solutions = system.sor_solve(omega=0.9)
    system.plot_convergence(residuals, "Convergence Plot")
    print("SOR Iterations:", len(residuals))
    print("SOR Solution:", x_sor)


if __name__ == "__main__":
    main()
