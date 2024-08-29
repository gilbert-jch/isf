function [d,val] = bf_optim(V)

%%
% [d] = bf_optim(V)
%
% Returns a direction d such that V'*d > 0 if any and d=[] otherwise.
% This is done by checking whether the optimal value of the linrar
% optimization problem
%
%   { min t
%   { V'*d+t*e >= 0
%   { t >= -1
%
% is negative (e is the vector of all ones).

% Version 1.1, September 1, 2023.
%
% HAL Id: hal-04124994 , version 1
% SWHID: swh:1:dir:0130daa41369602815de959df3fbe6de244c5763;origin=https://hal.archives-ouvertes.fr/hal-04124994;visit=swh:1:snp:1b642a2e1d9a9cf37b10b9c5487eccbad9b721e5;anchor=swh:1:rel:5208d7509416bc90534682e3d84f2adc08d183ad;path=/
%
% If you found this piece of software useful, please cite the paper:
% J.-P. Dussault, J.Ch. Gilbert, B. Plaquevent-Jourdain, 'On the
% B-differential of the componentwise minimum of two affine vectorial
% functions', 2023.

%-----------------------------------------------------------------------
%
% Authors:
% - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
% - Jean Charles Gilbert (INRIA, France),
% - Baptiste Plaquevent-Jourdain (INRIA & Univ. of Sherbrooke, Canada).
%
% Copyright 2023, INRIA (France) and Univ. of Sherbrooke (Canada).
%
% ISF is distributed under the terms of the Q Public License version
% 1.0.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the Q Public
% License version 1.0 for more details.
%
% You should have received a copy of the Q Public License version 1.0
% along with this program.  If not, see
% <https://doc.qt.io/archives/3.3/license.html>.
%
%-----------------------------------------------------------------------

% Get dimensions

  [n,p] = size(V);

% Compute lm

  options.Algorithm           = 'dual-simplex';
  options.Display             = 'off';	% final off iter
  options.Diagnostics         = 'off';
  options.OptimalityTolerance = 1.e-10;

  c       = zeros(n+1,1);
  c(n+1)  = 1;
  A       = -[[V';zeros(1,n)] ones(p+1,1)];
  b       = zeros(p+1,1);
  b(p+1)  = 1;

  [dt,val,flag] = linprog(c,A,b,[],[],[],[],[],options);
  d = dt(1:n);

  if flag ~= 1	% linprog failed
    d = [];
    return
  end

  return
