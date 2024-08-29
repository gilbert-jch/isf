function [info] = isf_memo_sc(p,qp,bvec,perm,info,options)

%%
% [info] = isf_memo_sc(p,qp,bvec,info)
%
% Recursive procedure to collect the infeasible s's by completing the
% infeasible bvec with 0 and 1's; qp = length(bvec), p is the desired
% length of the final bvec.

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
% <http://doc.trolltech.com/3.0/license.html>.
%
%-----------------------------------------------------------------------

% bvec has the right lenght p

  if qp == p
    info.nsc = info.nsc+1;
    if info.nsc > info.nscmax
      info.sc = [info.sc;zeros(info.nscmax,p)];
      info.nscmax = info.nscmax*2;
    end
    info.sc(info.nsc,perm) = bvec;
    return
  end

% Otherwise bvec is completed with 0 and 1's

  bvec = [bvec;0];
  [info] = isf_memo_sc(p,qp+1,bvec,perm,info,options);
  bvec(qp+1) = 1;
  [info] = isf_memo_sc(p,qp+1,bvec,perm,info,options);

  return
