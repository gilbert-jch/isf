function [W,colsel,colout,info] = isf_noncolin(V,info,values)

%%
% [W,colsel,colout,info] = isf_noncolin(V,info,values)
%
% This function selects in W the columns of V that are not colinear to
% another one. On return,
% - colsel gives the index of the selected columns of V,
% - colout(j), for j=setdiff(1:p,colsel), gives the index of the
%   colinear column of V.

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

% Dimension

  p = size(V,2);

% Small number

  eps2 = 100*eps;

% Initialization

  info.flag = values.success;

% Detect the columns of V that are not colinear to another one

  colsel = zeros(p,1);			% list of the selected noncolinear columns of V (at most p)
  colout = 1:p;				% pointers to colsel telling where are the colinear columns of V in the future W (positive if same sense)

  js = 1;				% pointer to colsel
  colsel(1) = 1;			% the first column of V is selected
  for jo = 2:p				% is the column V(:,jo) selected?
    v      = V(:,jo);			% compare with the selected columns of V
    select = true;
    for j = 1:js			% check whether V(jo) is colinear to a selected column
      cj = colsel(j);
      if norm(v-V(:,cj),Inf) <= eps2
        colout(jo) = cj;
        select = false;
        break;
      elseif norm(v+V(:,cj),Inf) <= eps2
        colout(jo) = -cj;
        select = false;
        break;
      end
    end
    if select
      js = js+1;
      colsel(js) = jo;
    end
  end

  W = V(:,colsel(1:js));		% define the new matrix V, named W, from the 'js' noncolinear columns of V
  colsel =  colsel(1:js);

  return
