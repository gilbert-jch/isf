function [twodesc,ratio,maxratiosI,minratiosJ] = isf_ntype(V,s,d,v,T)

%%
% [twodesc,ratio,maxratiosI,minratiosJ] = isf_ntype(V,s,d,v,T)
%

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

% Precision parameter

  eps3 = 1000*eps;

% Compute scalar products

  VTd  = V(:,T)'*d;
  VTv  = V(:,T)'*v;
  sVTv = s.*VTv;
  vTd  = v'*d;
  vn2  = v'*v;

% Determine whether s has 2 descendents for sure

  ratio  = -vTd/vn2;
  ratios = -VTd./VTv;
  I = find(sVTv > 0);
  J = find(sVTv < 0);

  if isempty(I)
    maxratiosI = -Inf;
    testI = true;
  else
    maxratiosI = max(ratios(I));
    testI = (maxratiosI < ratio - eps3);
  end
  if isempty(J)
    minratiosJ = +Inf;
    testJ = true;
  else
    minratiosJ = min(ratios(J));
    testJ = (ratio + eps3 < minratiosJ);
  end

  twodesc = testI & testJ;

  return
