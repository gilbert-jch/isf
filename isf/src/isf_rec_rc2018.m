function [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values)

%%
% [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values)
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

% Precision parameter

  eps5 = 1.e5*eps;	% small positive number to detect nonzero elements in the computed null space vecctor

% Sign vector s

  s = 2*bvec-1;

%-------------------------------------------------------------------------------------------------------------------------------
% Determine whether s has 1 or 2 descendents (depends on the options.dvnear0)
%-------------------------------------------------------------------------------------------------------------------------------

  v    = V(:,nvp);
  vTd  = v'*d;

  VTd  = V(:,1:nv)'*d;
  VTv  = V(:,1:nv)'*v;
  sVTv = s.*VTv;

  ratios = -VTd./VTv;
  I = find(sVTv >  eps);
  J = find(sVTv < -eps);

  % to avoid difficulties linked to an empty set I or J

  if isempty(I)
    maxratiosI = -Inf;
  else
    maxratiosI = max(ratios(I));
  end
  if isempty(J)
    minratiosJ = +Inf;
  else
    minratiosJ = min(ratios(J));
  end

  if (abs(vTd) < eps5) && (maxratiosI < -eps5) && (eps5 < minratiosJ)

% vTd near zero ==> 2 signs without optimization

    info.nb_npl(nvp) = info.nb_npl(nvp)+2;

    bvec = [bvec;1];
    if isempty(J)
      t = 1;
    else
      t  = 0.5*minratiosJ;
    end
    dp = d + t*v;
    dp = dp/norm(dp);		% ensure unit norm for d
    if verb2
      if (verb2 >= 3) || (p == nvp)
        [info] = isf_print(V,bvec,perm,'2+',dp,info,options,values);
        if info.flag; return; end
      end
    end
    if p > nvp
      [info] = isf_rec_rc2018(V,bvec,perm,dp,info,options,values);
      if info.flag; return; end
    else
      info.ns = info.ns+1;
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end

    bvec(nvp) = 0;
    if isempty(I)
      t = -1;
    else
      t = 0.5*maxratiosI;
    end
    dm = d + t*v;
    dm = dm/norm(dm);		% ensure unit norm for d
    if verb2
      if (verb2 >= 3) || (p == nvp)
        [info] = isf_print(V,bvec,perm,'2-',dm,info,options,values);
        if info.flag; return; end
      end
    end
    if p > nvp
      [info] = isf_rec_rc2018(V,bvec,perm,dm,info,options,values);
    else
      info.ns = info.ns+1;
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end

% Positive vTd ==> 1 single sign or 2 if optimization detects it

  elseif vTd > 0

    bvec = [bvec;1];
    if (verb2 >= 3) || (p == nvp)
      [info] = isf_print(V,bvec,perm,'>0',d,info,options,values);
      if info.flag; return; end
    end
    if p > nvp
      [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values);
      if info.flag; return; end
    else
      info.ns = info.ns+1;
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end
    info.nb_npl(nvp) = info.nb_npl(nvp)+1;

    % solve a LO problem to see whether bvec = [bvec;0] is a feasible binary vector

    bvec(nvp) = 0;
    time = tic;
    [d,~,info] = isf_feas((2*bvec(1:nv)-1)'.*V(:,1:nv),-v,d,info,options,values);
    info.cput_lop = info.cput_lop+toc(time);
    if info.flag; return; end
    if ~isempty(d)	% the linear inequality system is feasible
      if verb2
        [info] = isf_print(V,bvec,perm,'lp',d,info,options,values);
        if info.flag; return; end
      end
      if p > nvp
        [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values);
        if info.flag; return; end
      else
        info.ns = info.ns+1;
        if options.s
          if info.ns > info.nsmax
            info.s     = [info.s;zeros(info.nsmax,p)];
            info.nsmax = info.nsmax*2;
          end
          info.s(info.ns,perm) = bvec';
        end
      end
      info.nb_npl(nvp) = info.nb_npl(nvp)+1;
    elseif verb2 >= 3
      fprintf(fout2,'\n%0i | ',info.ns+1);
      [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
      fprintf(fout2,' | lp | infeasible system');
    end

% Negative vTd ==> 1 single sign or 2 if optimization detects it

  else

    bvec = [bvec;0];
    if (verb2 >= 3) || (p == nvp)
      [info] = isf_print(V,bvec,perm,'<0',d,info,options,values);
      if info.flag
        fprintf(fout2,'\n  v''*d = %12.5e is probably too close to 0',vTd);
        return
      end
    end
    if p > nvp
      [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values);
      if info.flag; return; end
    else
      info.ns = info.ns+1;
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end
    info.nb_npl(nvp) = info.nb_npl(nvp)+1;

    % solve a LO problem to see whether bvec = [bvec;1] is a feasible binary vector

    bvec(nvp) = 1;
    time = tic;
    [d,~,info] = isf_feas((2*bvec(1:nv)-1)'.*V(:,1:nv),v,d,info,options,values);
    info.cput_lop = info.cput_lop+toc(time);
    if info.flag; return; end
    if ~isempty(d)	% the linear inequality system is feasible
      if verb2
        [info] = isf_print(V,bvec,perm,'lp',d,info,options,values);
        if info.flag; return; end
      end
      if p > nvp
        [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values);
        if info.flag; return; end
      else
        info.ns = info.ns+1;
        if options.s
          if info.ns > info.nsmax
            info.s     = [info.s;zeros(info.nsmax,p)];
            info.nsmax = info.nsmax*2;
          end
          info.s(info.ns,perm) = bvec';
        end
      end
      info.nb_npl(nvp) = info.nb_npl(nvp)+1;
    elseif verb2 >= 3
      fprintf(fout2,'\n%0i | ',info.ns+1);
      [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
      fprintf(fout2,' | lp | infeasible system');
    end

  end

% Return

  return
