function [info] = isf_print(V,bvec,perm,type,d,info,options,values)

%%
% [info] = isf_print(V,bvec,perm,type,d,info,options,values)
%
% Print the binary vector bvec and do the appropirate verification if
% required.
%
% On entry
% - bvec = binary form of the sign vectors of the nv=legth(bvec)
%     considered vectors (columns) in V, which are
%     V(:,perm(1:length(bvec)));
% - type = 2character string giving the type of nodes;

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

% Default output

  info.flag = values.success;

% Get options

  fout   = options.fout;
  fout2  = options.fout2;
  verb2  = options.verb2;

% Get dimensions

  p  = size(V,2);	% total number of vectors
  nv = length(bvec);	% number of vectors with which bvec is associated

% Printings

  % Print the order number of the binary vector

  if verb2 >= 2
    if nv == p
      fprintf(fout2,'\n%0i | ',info.ns);
    else
      if info.ns == 0
        fprintf(fout2,'\n  | ');
      else
        fprintf(fout2,'\n%s | ',repmat(' ',1,floor(log10(info.ns))+1));
      end
    end
  end

  % Print the binary vector 'bvec' (and its complement if verb2 == 1) taking into account the permutation 'perm' of it indices

  [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
  
  % When verb >= 2, print also the calculus type and the direction

  if verb2 >= 2
    fprintf(fout2,' | %2s',type);
    fprintf(fout2,' |');
    if options.withd; fprintf(fout2,' %12.5e',d); end
  end

  % When verb >= 4, proceed with a verification of the direction d

  if options.withd && (verb2 >= 4)
    s = 2*bvec-1;
    sVd = zeros(p,1);
    if options.rc2018
      for i = 1:nv	% a loop since the compact multiplication is done columnwise and may destroy the sign of near zero values of vTd
        sVd(i) = s(i)*(V(:,perm(i))'*d);
      end
    else
      sVd = s.*(V(:,perm(1:nv))'*d);
    end
    if ~all(sVd > 0)
      info.flag = values.fail_on_sd_not_compatible;
      if fout == fout2
        fprintf(fout,' | failure\n\n  s.*(V(:,perm(1:%0i))''*d)\n  ^^^^^^^^^^^^^^^^^^^^^^^',nv);
        fprintf(fout,'\n  %23.16e',sVd);
        fprintf(fout,'\n');
      end
      return
    else
      if fout == fout2; fprintf(fout,' | success'); end
    end
  end

% Return

  return
