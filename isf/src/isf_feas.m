function [d,lm,info] = isf_feas(VS,sv,d0,info,options,values)

%%
% [d,info] = isf_feas(VS,sv,d0,info,options,values)
%
% Returns a nonempty d if and only if the system of strict linear
% inequalities
%
%   { VS'*d > 0
%   { sv'*d > 0
%
% is feasible for d (the returned d). In this system, VS is supposed to
% be a real nxp matrix and sv to be a real nx1 vector. This property is
% checked by using Linprog.

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

% Paramter

  eps2 = 100*eps;	% small number

% Set output

  d  = [];		% failure by default
  lm = [];

% Get dimensions

  [n,km] = size(VS);
  k      = km + 1;

% Compute d if any

  options_lp.Algorithm           = 'dual-simplex';	% 'dual-simplex' (useful to get many zeros in the dual variables) 'interior-point' (may claim that the problem is unbounded)
  options_lp.Display             = 'off';	% final off iter
  options_lp.Diagnostics         = 'off';
  options_lp.OptimalityTolerance = 1.e-9;
  options_lp.MaxIter             = 5000;

  % data of the LO problem
  %
  % min t for (d,t) in Rn x R
  % st  VS'*d >=  1               (-t in the RC paper)
  %     sv'*d >= -t
  %     t     >= -1

  c      = zeros(n+1,1);
  c(n+1) = 1;
  A      = -[[VS';sv';zeros(1,n)] [zeros(km,1);1;1]];
  b      = -ones(k+1,1);
  b(k)   = 0;
  b(k+1) = 1;

  % Compute a feasible starting point [d0;t0] from d0

  if options.rc2018
    minVSpd0 = Inf;
    for j = 1:km	% a loop since the compact multiplication is done columnwise and may destroy the sign of near zero values of vTd
      minVSpd0 = min(minVSpd0,VS(:,j)'*d0);
    end
  else
    minVSpd0 = min(VS'*d0);
  end
  if minVSpd0 <= 0
    info.flag = values.fail_on_technicality;
    return
  else
    d0 = (2/minVSpd0)*d0;	% ensures VS'*d0 >= 2
    t0 = max(-sv'*d0+1,0);	% ensures sv'*d >= -t0+1 and t0 >= 0
  end

% Ad0mb = A*[d0;t0]-b;

  % Run Linprog

  info.nb_losolve = info.nb_losolve+1;

  [dt,val,flag_lp,~,lambda] = linprog(c,A,b,[],[],[],[],[d0;t0],options_lp);

  % Diagnostics: system feasibility occurs only if val < 0

  if flag_lp ~= 1		% linprog failed
    if options.verb
      if flag_lp == 0
        fprintf('\n### Linprog: maximum number (%0i) of iterations reached',options_lp.MaxIter);
      elseif flag_lp == -2
        fprintf('\n### Linprog: no feasible point found');
      elseif flag_lp == -3
        fprintf('\n### Linprog: unbounded problem');
      elseif flag_lp == -4
        fprintf('\n### Linprog: NaN value encountered during execution of algorithm');
      elseif flag_lp == -5
        fprintf('\n### Linprog: both primal and dual problems are infeasible');
      elseif flag_lp == -7
        fprintf('\n### Linprog: magnitude of search direction became too small');
      else
        fprintf('\n### Linprog: outflag = %i',flag_lp);
      end
    end
    info.flag = values.fail_on_lp_solve;
    return
  elseif val < 0	% the system is feasible
    info.nb_feaslop = info.nb_feaslop+1;
    d = dt(1:n);
    d = d/norm(d);	% ensure unit norm for d
  else			% the system is infeasible
    info.nb_infeaslop = info.nb_infeaslop+1;
    if abs(lambda.ineqlin(k)-1) > eps2
      info.flag = values.fail_on_special_lp;
      return
    end
    lm = lambda.ineqlin(1:km);
  end

  info.flag = values.success;

  return
