function [info] = isf_sv_add(alpha,perm,info,options)

%%
% [info] = isf_sv_add(alpha,perm,info,options)
%
% If not already present, this function adds a stem vector in info.stem,
% using the null space vector alpha. One has V(:,perm(1:nv))*alpha == 0,
% where V is n x p matrix, nv = length(alpha), perm(1:nv) are the
% columns of V that accept alpha in its null space, V(:,perm(1:nv)) is
% of nullity one.

% To improve: it seems that s is used twice: once to define VS (V should
% suffice) and the second time below to define newsv.

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

% Dimensions

  p  = length(perm);	% total number of vectors in V
  nv = length(alpha);	% number of considered columns of V

% Parameter

  eps2 = 100*eps;	% small number

% Sign vector

  I          = find(abs(alpha)>eps2)';	% nonzero components of alpha (horizontal)
  s          = zeros(p,1);
  s(perm(I)) = sign(alpha(I));		% signed ordered indices in I of the linearly dependent vectors (horizontal)

% Get a restricted set of row of info.stems to scrutinize

  if (options.bestv > 0) || (ns == p)
    J = 1:info.nb_stems;
  else
    J = info.stem_zero_indices{p-nv};
  end

% See whether '±s' covers a stem vector in 'stems'

  adds = true;		% say whether 's' must be added to info.stems
  for j = J(:)'		% be sure that J(:)' is a row vector
    if abs(info.stems(j,:)*s) == info.stem_sizes(j)
      adds = false;
      break
    end
  end

% Add 's' to info.stems

  if adds

    % double the size of info.stems if it is filled (better than increasing it by one everytime)

    nsv = info.nb_stems;
    if nsv == info.max_stems
      info.stems      = [info.stems; zeros(info.max_stems,p)];
      info.stem_sizes = [info.stem_sizes; zeros(info.max_stems,1)];
      info.max_stems  = info.max_stems*2;
    end

    % add 's' to info.stems

    nsv                  = nsv+1;
    info.nb_stems        = nsv;
    info.nb_new_stems    = info.nb_new_stems+1;
    info.stems(nsv,:)    = s;
    info.stem_sizes(nsv) = sum(abs(s));

    % print if appropirate

    if options.verb2 >= 2

      fout2 = options.fout2;

      svec      = repmat(' ',p,1);	% binary vector as a string
      svec(s>0) = '1';
      svec(s<0) = '0';
      fprintf(fout2,'\n%s | %0s |      added in info.stems(%i)',repmat(' ',floor(log10(info.ns))+1,1),svec,nsv);

    end

  end

return
