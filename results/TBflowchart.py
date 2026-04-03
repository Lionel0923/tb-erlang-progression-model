import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import expon

# -----------------------------
# Parameters
# -----------------------------
rate = 2.0          # lambda
mean_dwell = 1/rate   # mean dwell time

# -----------------------------
# Theoretical PDF
# -----------------------------
x = np.linspace(0, 10, 500)
pdf = expon.pdf(x, scale=mean_dwell)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(7,5))

# exponential PDF (green)
plt.plot(x, pdf, color='green', linewidth=2, label='Simple(Exponential)')

# mean vertical line
plt.axvline(mean_dwell, color='black', linestyle='--', linewidth=2,
            label=f'Mean = {mean_dwell:.2f}')

plt.xlabel('Waiting Time')
plt.ylabel('Density')
plt.title('Simple(Exponential) Waiting Time Distribution')
plt.legend()

plt.tight_layout()
plt.show()