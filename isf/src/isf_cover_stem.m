function [I] = isf_cover_stem(bvec,perm,stems,stem_sizes,nb_stems,stem_zero_indices,bestv)

%%
% [I] = isf_cover_stem(bvec,perm,stems,stem_sizes,nb_stems,stem_zero_indices,bestv)
%
% Check whether the sign vector given by 'bvec' and 'perm' covers an
% element in the list 'stems' (if so, 'I' is the indices of those
% elements in 'stems', otherwise 'I' is empty). The list 'stems' is a
% matrix of size [nb_stems,p], each row being a signed vector with
% elements in {-1,0,+1} (0 meaning that the corresponding compoenent of
% the stem vector is not specified). One says that a sign vector 's'
% covers a stem vector 'stem' if stem(i) = s(i) for all i such that
% stem(i) ~= 0.

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

% Set dimensions

  p  = size(stems,2);
  ns = length(bvec);	% number of signs in the sign vector

% Get a restricted set of row of info.stems to scrutinize

  if (bestv > 0) || (ns == p)
    J = 1:nb_stems;
  else
    J = stem_zero_indices{p-ns};
  end

% Sign vector

  s             = zeros(p,1);
  s(perm(1:ns)) = 2*bvec-1;

% See whether '±s' covers a stem vector in 'stems'

% I = [];
% for j = J(:)'					% be sure that J(:)' is a row vector
%   if abs(stems(j,:)*s) == stem_sizes(j)	% then ±s covers stems(j,:)
%     I = j;
%     break
%   end
% end
% % This is less efficient than the following single instruction

  I = find(abs(stems(J,:)*s) == stem_sizes(J));

% Return

  return
