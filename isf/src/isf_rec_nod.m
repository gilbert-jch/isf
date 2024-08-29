function [info] = isf_rec_nod(V,bvec,perm,info,options,values)

%%
% [info] = isf_rec_nod(V,bvec,perm,info,options,values)
%
% Recursive procedure.

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

  fout2 = options.fout2;
  verb2 = options.verb2;

% Get dimensions

  p   = size(V,2);	% total number of vectors
  nv  = length(bvec);	% number of vectors with which bvec is associated
  nvp = nv+1;

%-------------------------------------------------------------------------------------------------------------------------------
% Try skp = +1 and -1 as the next sign
%-------------------------------------------------------------------------------------------------------------------------------

  bvec0 = bvec;
  isc   = [];		% force to check whether (s,+1) is feasible

  time = tic;
  if options.svsprod > 0
    svsprod0 = info.svsprod;	% memorize in a local variable the value of info.svsprod on entry (the latter is modified in this function and by the recursive calls)
  end
  info.cput_cover = info.cput_cover+toc(time);

  for skp = [1,0]

    bvec = [bvec0;skp];

    if skp == 1
      skp_string = '+1';
    else
      skp_string = '-1';
    end

    % when options.svsprod > 0, compute info.svsprod for the current bvec, whatever isc is

    time = tic;

    if options.svsprod > 0
      if skp == 1							% to avoid a multiplication vector*(±1)
        svsprod = svsprod0 + info.stems(1:info.nb_stems,perm(nvp));	% value of svsprod at the current S-tree level (to get it on return from the recursive process)
      else
        svsprod = svsprod0 - info.stems(1:info.nb_stems,perm(nvp));	% value of svsprod at the current S-tree level (to get it on return from the recursive process)
      end
      info.svsprod = svsprod;						% for the recursive call
    end

    % see whether bvec is feasible

    if isempty(isc)
      info.nb_svdetect = info.nb_svdetect+1;
      if options.svsprod == 0
        isc = isf_cover_stem(bvec,perm,info.stems,info.stem_sizes,info.nb_stems,info.stem_zero_indices,options.bestv);
      else
        isc = find(abs(svsprod) == info.stem_sizes);			% indices of the stem vectors covered by the current sign vector, if any
      end
    else		% (s,+1) is infeasible, hence (s,-1) is feasible
      isc = [];
    end

    info.cput_cover = info.cput_cover+toc(time);

    if ~isempty(isc)		% bvec is infeasible

      if verb2 >= 3
	if info.ns == 0
          fprintf(fout2,'\n  | ');
	else
          fprintf(fout2,'\n%s | ',repmat(' ',1,floor(log10(info.ns))+1));
	end
        [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
        fprintf(fout2,' | %s | infeasible system (covering info.stems(%i))',skp_string,isc(1));
      end

    else			% bvec is feasible

      if nvp == p; info.ns = info.ns+1; end
      info.nb_npl(nvp) = info.nb_npl(nvp)+1;
      if (verb2 >= 3) || (p == nvp)
        [info] = isf_print(V,bvec,perm,skp_string,[],info,options,values);
        if info.flag; return; end
      end

      % recursive call

      if p > nvp
        [info] = isf_rec_nod(V,bvec,perm,info,options,values);
        if info.flag; return; end
      else
        if options.s; info.s(info.ns,perm) = bvec'; end
      end

    end

  end

  return
