import numpy as np
import os
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt

current_dir = os.path.dirname(os.path.abspath(__file__))
picture_dir = os.path.join(current_dir, "..", "picture")
os.makedirs(picture_dir, exist_ok=True)


# Original data
x_data = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
y_original = np.array([3.5, 4.2, 3.8, 4.5, 4.0, 3.6])

# 1. Build original spline
cs_original = CubicSpline(x_data, y_original, bc_type="natural")

# 2. Local modification: add 0.2m in region x=30 to x=50
y_modified = y_original.copy()

# Determine nodes to modify (Method 1: directly modify relevant nodes)
modify_indices = [1, 2]  # Index 1(x=20) and 2(x=40) affect [30,50] region
y_modified[modify_indices] += 0.2  # Approximate adjustment

# More precise method: insert new nodes at 30 and 50 (Method 2: node refinement)
x_extended = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0])
y_extended = np.array([3.5, 4.2, 4.0, 4.5, 4.0, 3.6])

# 3. Build updated spline
cs_modified = CubicSpline(x_extended, y_extended, bc_type="natural")

# Check continuity conditions
x_test = np.linspace(0, 100, 1000)
y_orig = cs_original(x_test)
y_mod = cs_modified(x_test)

# Calculate first derivative continuity
dy_orig = cs_original(x_test, 1)
dy_mod = cs_modified(x_test, 1)

# Calculate second derivative continuity
d2y_orig = cs_original(x_test, 2)
d2y_mod = cs_modified(x_test, 2)


def analyze_impact(cs_orig, cs_mod, x_range):
    """Analyze the impact of modifications on the curve"""
    x_fine = np.linspace(x_range[0], x_range[1], 500)

    # Calculate changes
    delta_y = cs_mod(x_fine) - cs_orig(x_fine)
    delta_dy = cs_mod(x_fine, 1) - cs_orig(x_fine, 1)
    delta_d2y = cs_mod(x_fine, 2) - cs_orig(x_fine, 2)

    return x_fine, delta_y, delta_dy, delta_d2y


# Analyze impact in different regions
regions = {
    "Unaffected Upstream Region": [0, 25],
    "Modification Boundary Region": [25, 35],
    "Core Modification Region": [35, 55],
    "Affected Downstream Region": [55, 100],
}

plt.figure(figsize=(12, 10))

# 1. Comparison of original and modified curves
plt.subplot(3, 1, 1)
plt.plot(x_test, y_orig, "b-", label="Original Curve", linewidth=2)
plt.plot(x_test, y_mod, "r--", label="Modified Curve", linewidth=2)
plt.plot(x_data, y_original, "bo", markersize=8, label="Original Measurement Points")
plt.plot(x_extended, y_extended, "rx", markersize=6, label="Modified Points")
plt.axvspan(30, 50, alpha=0.2, color="yellow", label="Modified Region")
plt.ylabel("Draft Depth (m)")
plt.title("Cubic Spline Interpolation: Local Modification Comparison")
plt.legend()
plt.grid(True, alpha=0.3)

# 2. Draft depth change
plt.subplot(3, 1, 2)
delta_y = y_mod - y_orig
plt.plot(x_test, delta_y, "g-", linewidth=2)
plt.axvspan(30, 50, alpha=0.2, color="yellow")
plt.axhline(y=0, color="k", linestyle="--", alpha=0.5)
plt.ylabel("Draft Change (m)")
plt.title("Impact Propagation of Local Modification")
plt.grid(True, alpha=0.3)

# 3. Slope change
plt.subplot(3, 1, 3)
delta_dy = cs_modified(x_test, 1) - cs_original(x_test, 1)
plt.plot(x_test, delta_dy, "m-", linewidth=2)
plt.axvspan(30, 50, alpha=0.2, color="yellow")
plt.axhline(y=0, color="k", linestyle="--", alpha=0.5)
plt.xlabel("Position x (m)")
plt.ylabel("Slope Change")
plt.title("First Derivative Change")
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(
    os.path.join(picture_dir, "local_modification.png"), dpi=300, bbox_inches="tight"
)
plt.show()
