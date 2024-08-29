function [info] = isf_print_bvec(fout,verb,p,bvec,perm,info)

%%
% [info] = isf_print_bvec(fout,verb,p,bvec,perm,info)
%
% Print the binary vector 'bvec' on channel 'fout' , taking into account
% the permutation 'perm' of its indices.

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

% Print bvec

  nb = length(bvec);

% Otherwise, print bvec

  if (verb >= 2) || (nb == p)
    svec = repmat(' ',p,1);		% binary vector as a character array
    svec(perm(bvec==1))='1';
    svec(perm(bvec==0))='0';
    if verb >= 2; fprintf(fout,'%0s',svec); end
  end

% When verb == 1, print also the binary complementary of bvec

  if (verb == 1) && (nb == p)
    cvec = ones(p,1)-bvec;		% complementary binary vector
    svec = repmat(' ',p,1);		% binary vector as a character array
    svec(perm(cvec==1))='1';
    svec(perm(cvec==0))='0';
    fprintf(fout,'\n  %0s',svec);
  end

% Print bvec

  return
