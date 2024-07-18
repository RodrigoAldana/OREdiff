# Copyright (C) 2022 Rodrigo Aldana-López
# <rodrigo.aldana.lopez at {gmail,intel} dot com> (Intel Corporation)
# For more information see <https://github.com/RodrigoAldana/OREdiff/blob/main/README.md>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import matplotlib.pyplot as plt

class Differentiator:
    def __init__(self, kbar):
        self.buffer = np.zeros(kbar + 1)
        self.kbar = kbar

    def update_diff(self, f, k, Delta, L):
        self.buffer[:-1] = self.buffer[1:]
        self.buffer[-1] = f
        Nhat = 0

        for l in range(2, min(k, self.kbar) + 1):
            for j in range(1, l + 1):
                Q = self.buffer[-j] - self.buffer[-1] + (self.buffer[-1] - self.buffer[-l]) * (j / l)
                Nlj = (abs(Q) - L * Delta**2 * j * (l - j) / 2) / 2
                if Nlj > Nhat:
                    Nhat = Nlj

        if Nhat == 0:
            lhat = 1
        else:
            lhat = min(k, self.kbar, int(np.ceil(2 * np.sqrt(Nhat / L) / Delta)))

        y = (self.buffer[-1] - self.buffer[-lhat]) / (lhat * Delta)
        return y, Nhat

def update_filter(y, input, Delta, gamma=1.96):
    if abs(input - y) <= gamma * Delta:
        sat = input - y
    else:
        sat = gamma * Delta * np.sign(input - y)
    ynew = y + sat
    return ynew

if __name__ == '__main__':
    np.random.seed(0)

    Delta = 0.01
    L = 1
    R = 1
    Nbar = 0.08
    kbar = int(np.ceil(np.sqrt(2 * Nbar / (L * Delta**2)) + 1))

    T = 20
    t = np.arange(0, T + Delta, Delta)

    c = 0.5
    f = L * np.sin(t)
    df = L * np.cos(t)

    noise = 2 * Nbar * (np.random.rand(len(t)) - 0.5)
    y_diff = np.zeros(len(t))
    y_filter = np.zeros(len(t))
    lhat = np.zeros(len(t))
    Nhat = np.zeros(len(t))

    differentiator = Differentiator(kbar)

    for k in range(1, len(t)):
        if k % 100 == 0:
            print(f"{k / len(t):.2%}")
        input_signal = f[k] + noise[k]
        y_diff[k], Nhat[k] = differentiator.update_diff(input_signal, k, Delta, L)
        y_filter[k] = update_filter(y_filter[k-1], y_diff[k], Delta)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    ax1.plot(t, df, 'b--', linewidth=2, label='derivative of f')
    ax1.plot(t, y_diff, color=[0.5, 0.5, 0.5], linewidth=1, label='optimal differentiator')
    ax1.plot(t, y_filter, 'r', linewidth=2, label='filtered differentiator')
    ax1.set_xlabel('time (seconds)')
    ax1.set_ylabel('derivative, estimation')
    ax1.legend()
    ax1.grid(True)

    A = 2 * np.sqrt(2 * Nbar * L) + L * Delta / 2
    ax2.plot(t, np.abs(y_diff - df), color=[0.5, 0.5, 0.5], linewidth=1, label='Optimal differentiator')
    ax2.plot(t, np.abs(y_filter - df), 'r', linewidth=2, label='Filtered differentiator')
    ax2.plot([0, t[-1]], [A, A], 'k--', linewidth=2, label='Worst-case accuracy bound')
    ax2.set_xlabel('time (seconds)')
    ax2.set_ylabel('differentiation error')
    ax2.legend()
    ax2.axis([0, T, 0, A * 1.1])
    ax2.grid(True)

    plt.tight_layout()
    plt.show()

