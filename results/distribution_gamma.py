import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gamma

# -----------------------------
# Parameters
# -----------------------------
mean_dwell = 0.5   # mean dwell time
k1 = 2
k2 = 50

# rates adjusted so mean is the same
rate1 = k1 / mean_dwell
rate2 = k2 / mean_dwell

# x range
x = np.linspace(0, 6, 500)

# Erlang PDFs (gamma distribution with integer shape)
pdf_k2 = gamma.pdf(x, a=k1, scale=1/rate1)
pdf_k50 = gamma.pdf(x, a=k2, scale=1/rate2)

# -----------------------------
# Plot
# -----------------------------
plt.figure(figsize=(7,5))

plt.plot(x, pdf_k2, linewidth=2, label='Gamma-distributed k=2')
plt.plot(x, pdf_k50, linewidth=2, label='Gamma-distributed k=50')

# mean vertical line
plt.axvline(mean_dwell, color='black', linestyle='--',
            label=f'Mean = {mean_dwell}')

plt.xlabel('Waiting Time')
plt.ylabel('Density')
plt.title('Gamma Distributed Dwell Time Distributions')
plt.legend()

plt.tight_layout()
plt.show()