function [y, Nhat] = update_diff(f, k, kbar, Delta, L)
%% update_diff - Differentiator for a signal with noise amplitude estimation
% Inputs:
%    f     - Signal to be differentiated
%    k     - Time instant
%    kbar  - Window length for the memory
%    Delta - Discretization step
%    L     - Lipschitz constant for the derivative of y
%
% Outputs:
%    y     - Differentiated signal
%    Nhat  - Noise amplitude estimation
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
    persistent buffer;
    if isempty(buffer)
        buffer = zeros(1, kbar + 1);
    end

    buffer(1:end-1) = buffer(2:end);
    buffer(end) = f;
    Nhat = 0;

    for l = 2:min(k, kbar)
        for j = 1:l
            Q = buffer(end-j) - buffer(end) + (buffer(end) - buffer(end-l)) * (j / l);
            Nlj = (abs(Q) - L * Delta^2 * j * (l - j) / 2) / 2;
            if Nlj > Nhat
                Nhat = Nlj;
            end
        end
    end

    if Nhat == 0
        lhat = 1;
    else
        lhat = min([k, kbar, ceil(2 * sqrt(Nhat / L) / Delta)]);
    end

    y = (buffer(end) - buffer(end - lhat)) / (lhat * Delta);
end
