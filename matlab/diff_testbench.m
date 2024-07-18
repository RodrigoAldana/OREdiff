%%
% Copyright (C) 2022 Rodrigo Aldana-López
% <rodrigo.aldana.lopez at {gmail,intel} dot com> (Intel Corporation)
% For more information see <https://github.com/RodrigoAldana/OREdiff/blob/main/README.md>
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%
clc;
clear all;
close all;

rng('default');

Delta = 0.01;
L = 1;
R = 1;
Nbar = 0.08;
kbar = ceil(sqrt(2 * Nbar / (L * Delta^2)) + 1);

T = 20;
t = 0:Delta:T;

c = 0.5;
f = L * sin(t);
df = L * cos(t);

noise = 2 * Nbar * (rand(size(t)) - 0.5);
y_diff = zeros(size(t));
y_filter = zeros(size(t));
lhat = zeros(size(t));
Nhat = zeros(size(t));

for k = 2:numel(t)
    if mod(k, 100) == 0
        disp(k / numel(t));
    end
    input = f(k) + noise(k);
    [y_diff(k), Nhat(k)] = update_diff(input, k, kbar, Delta, L);
    y_filter(k) = update_filter(y_filter(k-1), y_diff(k), Delta);
end

fig = figure();
subplot(2,1,1)

plot(t, df, 'b--','LineWidth', 2); hold on;
plot(t, y_diff, 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1); hold on;
plot(t, y_filter, 'r', 'LineWidth', 2);
xlabel('time (seconds)')
ylabel('derivative, estimation')
legend('derivative of f', 'optimal differentiator', 'filtered differentiator')
grid on;

subplot(2,1,2)
A = 2 * sqrt(2 * Nbar * L) + L * Delta / 2;
plot(t, abs(y_diff - df), 'Color', [0.5, 0.5, 0.5], 'LineWidth', 1); hold on;
plot(t, abs(y_filter - df), 'r', 'LineWidth', 2);
plot([0, t(end)], [A, A], 'k--', 'LineWidth', 2);
xlabel('time (seconds)')
ylabel('differentiation error')
legend('Optimal differentiator', 'Filtered differentiator', 'Worst-case accuracy bound')
axis([0, T, 0, A * 1.1]);
grid on;

fig.Renderer = 'Painters';
