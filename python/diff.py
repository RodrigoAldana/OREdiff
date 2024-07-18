import numpy as np
import matplotlib.pyplot as plt

class Differentiator:
    def __init__(self, L, Delta, Nbar, gamma=1.96):
        self.L = L
        self.Delta = Delta
        self.Nbar = Nbar
        self.gamma = gamma
        self.kbar = int(np.ceil(np.sqrt(2 * Nbar / (L * Delta**2)) + 1))
        self.buffer = np.zeros(self.kbar + 1)
        self.y_filter = 0
        self.k = 0

    def diff(self, f, use_filter=False):
        self.buffer[:-1] = self.buffer[1:]
        self.buffer[-1] = f
        self.k += 1
        Nhat = 0

        for l in range(2, min(self.k, self.kbar) + 1):
            for j in range(1, l + 1):
                Q = self.buffer[-j] - self.buffer[-1] + (self.buffer[-1] - self.buffer[-l]) * (j / l)
                Nlj = (abs(Q) - self.L * self.Delta**2 * j * (l - j) / 2) / 2
                if Nlj > Nhat:
                    Nhat = Nlj

        if Nhat == 0:
            lhat = 1
        else:
            lhat = min(self.k, self.kbar, int(np.ceil(2 * np.sqrt(Nhat / self.L) / self.Delta)))

        y = (self.buffer[-1] - self.buffer[-lhat]) / (lhat * self.Delta)
        
        if use_filter:
            y = self.filter(y)
        
        return y

    def filter(self, y):
        if abs(y - self.y_filter) <= self.gamma * self.Delta:
            sat = y - self.y_filter
        else:
            sat = self.gamma * self.Delta * np.sign(y - self.y_filter)
        self.y_filter += sat
        return self.y_filter

if __name__ == '__main__':
    np.random.seed(0)

    Delta = 0.01
    L = 1
    Nbar = 0.08
    gamma = 1.96

    T = 20
    t = np.arange(0, T + Delta, Delta)

    f = L * np.sin(t)
    df = L * np.cos(t)

    noise = 2 * Nbar * (np.random.rand(len(t)) - 0.5)
    y_diff = np.zeros(len(t))
    y_diff_filtered = np.zeros(len(t))

    differentiator = Differentiator(L, Delta, Nbar, gamma)
    differentiatorF = Differentiator(L, Delta, Nbar, gamma)

    for k in range(1, len(t)):
        if k % 100 == 0:
            print(f"{k / len(t):.2%}")
        input_signal = f[k] + noise[k]
        y_diff[k] = differentiator.diff(input_signal, use_filter=False)
        y_diff_filtered[k] = differentiatorF.diff(input_signal, use_filter=True)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    ax1.plot(t, df, 'b--', linewidth=2, label='derivative of f')
    ax1.plot(t, y_diff, color=[0.5, 0.5, 0.5], linewidth=1, label='optimal differentiator')
    ax1.plot(t, y_diff_filtered, 'r', linewidth=2, label='filtered differentiator')
    ax1.set_xlabel('time (seconds)')
    ax1.set_ylabel('derivative, estimation')
    ax1.legend()
    ax1.grid(True)

    A = 2 * np.sqrt(2 * Nbar * L) + L * Delta / 2
    ax2.plot(t, np.abs(y_diff - df), color=[0.5, 0.5, 0.5], linewidth=1, label='Optimal differentiator')
    ax2.plot(t, np.abs(y_diff_filtered - df), 'r', linewidth=2, label='Filtered differentiator')
    ax2.plot([0, t[-1]], [A, A], 'k--', linewidth=2, label='Worst-case accuracy bound')
    ax2.set_xlabel('time (seconds)')
    ax2.set_ylabel('differentiation error')
    ax2.legend()
    ax2.axis([0, T, 0, A * 1.1])
    ax2.grid(True)

    plt.tight_layout()
    plt.show()
