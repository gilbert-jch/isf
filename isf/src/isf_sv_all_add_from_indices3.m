function [info,newstem] = isf_sv_all_add_from_indices3(V,cols,info)

%%
% [info,newstem] = isf_sv_all_add_from_indices3(V,cols,info)
%
% From a circuit of V, whose columns are included in those given in
% 'cols', construct a sign-vector (a list of signed integer, whose
% absolute value is in increasing) and put it in info.stems is
% appropriate. Appropriate means that either it is not in info.stems or
% it replaces another circuit in info.stems.

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

% Precision parameter

  eps5 = 1.e5*eps;	% small positive number to detect nonzero elements in the computed null space vecctor (100*eps is too small on 'data_ratio_20_5_7.txt')

% See whether the null space has dimension 1

  Z = null(V(:,cols));			% it is a matrix, whose columns form a basis of the kernel of V(:,cols) (vertical columns)
  if size(Z,2) ~= 1			% the nullity is ~= 1
    newstem = false;
    return
  end
  newstem = true;			% a stem vector is found

  I             = find(abs(Z)>eps5);	% nonzero components of the null space
  stem          = zeros(1,p);
  stem(cols(I)) = sign(Z(I));		% one of the two stem vectors
  if Z(I(1)) < 0
    stem = -stem;			% force the first nonzero component to be +1 (useful for identifying duplicates with 'unique')
  end

% Double the size of info.stems if it is filled (better than increasing it by one everytime)

  if info.nb_stems == info.max_stems
    info.stems      = [info.stems; zeros(info.max_stems,p)];
    info.max_stems  = info.max_stems*2;
  end

% Memorize stem in info.stems

  info.nb_stems = info.nb_stems+1;

  info.stems(info.nb_stems,:) = stem;

  return
