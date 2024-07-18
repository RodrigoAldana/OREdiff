function ynew = update_filter(y, input, Delta, gamma)
%% update_filter - Filters a signal using a specified Lipschitz constant
% Inputs:
%    y     - Previous value of the output
%    input - Signal to be filtered
%    Delta - Time step
%    gamma - Lipschitz constant for the output
%
% Outputs:
%    ynew  - Filtered output signal
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
    % Set the default value of gamma if not provided
    if nargin < 4
        gamma = 1.96;
    end

    if abs(input - y) <= gamma * Delta
        sat = input - y;
    else
        sat = gamma * Delta * sign(input - y);
    end
    ynew = y + sat;
end
