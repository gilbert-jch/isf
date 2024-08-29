function [b,p] = isf_btree(b,p)

%%
% [b] = isf_btree(b)
%
% Gives the next element in a binary tree from the current one given by
% b(1:p): b is a binary vector whose length is the depth of the tree and
% p is a pointer in [1:length(b)] giving the level of the tree of the
% current position. The ones come before the zeros. The tree is
% completely explored when isf_btree returns p=0. The first node down
% the root is described by b = [1 0 ... 0] and p = 1.

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

  pmax = length(b);

  while true	% relax; it is safe

    % if possible, go down with a 1 and return

    if p < pmax
      p    = p+1;
      b(p) = 1;
      return
    end

    % here p == length(b); hence backtrack to find the first 1, replace it by 0 and loop

    while (p > 0) && (b(p) == 0); p = p-1; end	% backtrack to find the first 1
    if p == 0; return; end			% end of the tree exploration
    b(p) = 0;					% replace 1 by 0 and loop

  end
