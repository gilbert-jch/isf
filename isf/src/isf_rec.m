function [info] = isf_rec(V,bvec,perm,d,info,options,values)

%%
% [info] = isf_rec(V,bvec,perm,d,info,options,values)
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

  eps3 = 1.e3*eps;	% small positive number to detect nonzero elements in the computed null space vecctor

% Sign vector s

  s = 2*bvec-1;

%-------------------------------------------------------------------------------------------------------------------------------
% Determine the next v; this operation must be done for each {s_i v_i: i in T}, since although the signs s_i play no role (with
% the adopted technique), d may change (and this has an impoct on the choice of the new vector with the adopted techniques)
%-------------------------------------------------------------------------------------------------------------------------------

  T = perm(1:nv);		% the index set of the vectors v used so far

  % determine the index 'ivp' of the next v (recall that the columns of V have unit norm)

  if options.bestv == 0		% choice (C0)

    ivp = perm(nvp);

  else

    Tc = setdiff(1:p,T);		% complementary set of T (horizontal)

    if options.bestv == 1	% choice (C1)

      [~,I] = max(abs(V(:,Tc)'*d));
      ivp   = Tc(I(1));			% index of the selected vector

    elseif options.bestv == 2	% choice (C2)

      dd    = sum(V(:,T),2);		% sum of the vectors considered so far
      [~,I] = max(abs(V(:,Tc)'*dd));
      ivp   = Tc(I(1));			% index of the selected vector

    elseif options.bestv == 3	% choice (C3)

      Tb = [];				% list of indices of vectors that may have a single descendant
      for i = Tc
        v = V(:,i);
        [twodesc,ratio,maxratiosI,minratiosJ] = isf_ntype(V,s,d,v,T);
        if ~twodesc
          Tb = [Tb,i];
        end
      end
      if isempty(Tb); Tb = Tc; end
      [~,I] = max(abs(V(:,Tb)'*d));
      ivp   = Tb(I(1));			% index of the selected vector

    end

    % index of the selected vector

    perm(nvp) = ivp;			% memorize ivp in perm

  end

  v   = V(:,ivp);
  vTd = v'*d;

%-------------------------------------------------------------------------------------------------------------------------------
% Determine whether s has 1 or 2 descendents (depends on the options.dvnear0)
%-------------------------------------------------------------------------------------------------------------------------------

  if options.dvnear0 || (abs(vTd) < eps3)
    [twodesc,ratio,maxratiosI,minratiosJ] = isf_ntype(V,s,d,v,T);
  else
    twodesc = false;
  end

%-------------------------------------------------------------------------------------------------------------------------------
% vTd near zero ==> s has 2 descendants (without having to use an optimization problem to decide)
%-------------------------------------------------------------------------------------------------------------------------------

  if twodesc

    % descendant (s,+1)

    if nvp == p; info.ns = info.ns+1; end	% can one have nvp == p here?
    info.nb_npl(nvp) = info.nb_npl(nvp)+1;
    bvec = [bvec;1];
    if minratiosJ == +Inf
      t = ratio+1;
    else
      t  = 0.5*(ratio+minratiosJ);
    end
    dp = d + t*v;
    dp = dp/norm(dp);		% ensure unit norm for d
    if verb2
      if (verb2 >= 3) || (p == nvp)
        [info] = isf_print(V,bvec,perm,'2+',dp,info,options,values);
        if info.flag; return; end
      end
    end
    if nvp < p
      % compute info.svsprod for the current S-tree node (described by bvec)
      time = tic;
      if options.svsprod > 0
        info.svsprod = svsprod0 + info.stems(1:info.nb_stems,perm(nvp));
      end
      info.cput_cover = info.cput_cover+toc(time);
      % recursive call
      [info] = isf_rec(V,bvec,perm,dp,info,options,values);
      if info.flag; return; end
    else
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end

    % descendant (s,-1)

    if nvp == p; info.ns = info.ns+1; end
    info.nb_npl(nvp) = info.nb_npl(nvp)+1;
    bvec(nvp) = 0;
    if maxratiosI == -Inf
      t = ratio-1;
    else
      t = 0.5*(ratio+maxratiosI);
    end
    dm = d + t*v;
    dm = dm/norm(dm);		% ensure unit norm for d
    if verb2
      if (verb2 >= 3) || (p == nvp)
        [info] = isf_print(V,bvec,perm,'2-',dm,info,options,values);
        if info.flag; return; end
      end
    end
    if nvp < p
      % compute info.svsprod for the current S-tree node (described by bvec)
      time = tic;
      if options.svsprod > 0
        info.svsprod = svsprod0 - info.stems(1:info.nb_stems,perm(nvp));
      end
      info.cput_cover = info.cput_cover+toc(time);
      % recursive call
      [info] = isf_rec(V,bvec,perm,dm,info,options,values);
    else
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end

%-------------------------------------------------------------------------------------------------------------------------------
% Positive vTd ==> 1 single sign vector or 2 if optimization detects it
%-------------------------------------------------------------------------------------------------------------------------------

  elseif vTd > 0

    % one sign is trivially feasible

    if nvp == p; info.ns = info.ns+1; end
    bvec = [bvec;1];
    if (verb2 >= 3) || (p == nvp)
      [info] = isf_print(V,bvec,perm,'>0',d,info,options,values);
      if info.flag; return; end
    end

    if nvp < p
      % compute info.svsprod for the current S-tree node (described by bvec)
      time = tic;
      if options.svsprod > 0
        info.svsprod = svsprod0 + info.stems(1:info.nb_stems,perm(nvp));
      end
      info.cput_cover = info.cput_cover+toc(time);
      % recursive call
      [info] = isf_rec(V,bvec,perm,d,info,options,values);
      if info.flag; return; end
    else
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end
    info.nb_npl(nvp) = info.nb_npl(nvp)+1;

    % consider now the other sign

    bvec(nvp) = 0;

    % compute info.svsprod for the current S-tree node (described by bvec)

    time = tic;
    if options.svsprod > 0
      info.svsprod = svsprod0 - info.stems(1:info.nb_stems,perm(nvp));
    end
    info.cput_cover = info.cput_cover+toc(time);

    % if option.sv > 0, verify whether ±s covers an element in info.stems, in which case s is infeasible

    if options.sv == 0
      isc = [];		% this means that ±s does not cover any stem vector (actually there is no stem vector when options.sv == 0)
    else
      time = tic;
      isc = isf_cover_stem(bvec,perm,info.stems,info.stem_sizes,info.nb_stems,info.stem_zero_indices,options.bestv);
      info.cput_cover = info.cput_cover+toc(time);
    end

    if ~isempty(isc)	% ±s covers the stem vector 'isc', hence bvec is an infeasible binary vector

      info.nb_svdetect = info.nb_svdetect+1;
      if verb2 >= 3
        fprintf(fout2,'\n%s | ',repmat(' ',1,floor(log10(info.ns))+1));
        [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
        fprintf(fout2,' | sv | infeasible system (covering info.stems(%i))',isc);
      end
      if options.sc; [info] = isf_memo_sc(p,nvp,bvec,perm,info,options); end

    else		% solve a LO problem to see whether bvec = [bvec;0] is a feasible binary vector

      time = tic;
      [d,lm,info] = isf_feas((2*bvec(1:nv)-1)'.*V(:,T),-v,d,info,options,values);
      info.cput_lop = info.cput_lop+toc(time);
      if info.flag; return; end

      if ~isempty(d)	% the linear inequality system is feasible
        if nvp == p; info.ns = info.ns+1; end
        if verb2
          [info] = isf_print(V,bvec,perm,'lp',d,info,options,values);
          if info.flag; return; end
        end
        if p > nvp
          [info] = isf_rec(V,bvec,perm,d,info,options,values);
          if info.flag; return; end
        else
          if options.s
            if info.ns > info.nsmax
              info.s     = [info.s;zeros(info.nsmax,p)];
              info.nsmax = info.nsmax*2;
            end
            info.s(info.ns,perm) = bvec';
          end
        end
        info.nb_npl(nvp) = info.nb_npl(nvp)+1;
      else		% the linear inequality system is infeasible
        if verb2 >= 3
          fprintf(fout2,'\n%s | ',repmat(' ',1,floor(log10(info.ns))+1));
          [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
          fprintf(fout2,' | lp | infeasible system');
        end
        if options.sc; [info] = isf_memo_sc(p,nvp,bvec,perm,info,options); end
        if options.sv >= 2 	% add a new stem vector to info.stems
          [info] = isf_sv_add((2*bvec-1).*[lm;1],perm,info,options);	% update info.stems by adding a sign vector if appropriate
        end
      end

    end

%-------------------------------------------------------------------------------------------------------------------------------
% Negative vTd ==> 1 single sign vector or 2 if optimization detects it
%-------------------------------------------------------------------------------------------------------------------------------

  else

    % one sign is trivially feasible

    if nvp == p; info.ns = info.ns+1; end
    bvec = [bvec;0];
    if (verb2 >= 3) || (p == nvp)
      [info] = isf_print(V,bvec,perm,'<0',d,info,options,values);
      if info.flag; return; end
    end
    if nvp < p
      % compute info.svsprod for the current S-tree node (described by bvec)
      time = tic;
      if options.svsprod > 0
        info.svsprod = svsprod0 - info.stems(1:info.nb_stems,perm(nvp));
      end
      info.cput_cover = info.cput_cover+toc(time);
      % recursive call
      [info] = isf_rec(V,bvec,perm,d,info,options,values);
      if info.flag; return; end
    else
      if options.s
        if info.ns > info.nsmax
          info.s     = [info.s;zeros(info.nsmax,p)];
          info.nsmax = info.nsmax*2;
        end
        info.s(info.ns,perm) = bvec';
      end
    end
    info.nb_npl(nvp) = info.nb_npl(nvp)+1;

    % consider now the other sign

    bvec(nvp) = 1;

    % compute info.svsprod for the current S-tree node (described by bvec)

    time = tic;
    if options.svsprod > 0
      info.svsprod = svsprod0 + info.stems(1:info.nb_stems,perm(nvp));
    end
    info.cput_cover = info.cput_cover+toc(time);

    % if option.sv > 0, verify whether ±s covers an element in info.stems, in which case s is infeasible

    if options.sv == 0
      isc = [];		% this means that ±s does not cover any stem vector (actually there is no stem vector when options.sv == 0)
    else
      time = tic;
      isc = isf_cover_stem(bvec,perm,info.stems,info.stem_sizes,info.nb_stems,info.stem_zero_indices,options.bestv);
      info.cput_cover = info.cput_cover+toc(time);
    end

    if ~isempty(isc)	% ±s covers the stem vector 'isc', hence bvec is an infeasible binary vector

      info.nb_svdetect = info.nb_svdetect+1;
      if verb2 >= 3
        fprintf(fout2,'\n%s | ',repmat(' ',1,floor(log10(info.ns))+1));
        [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
        fprintf(fout2,' | sv | infeasible system (covering info.stems(%i))',isc);
      end
      if options.sc; [info] = isf_memo_sc(p,nvp,bvec,perm,info,options); end

    else		% solve a LO problem to see whether bvec = [bvec;1] is a feasible binary vector

      time = tic;
      [d,lm,info] = isf_feas((2*bvec(1:nv)-1)'.*V(:,T),v,d,info,options,values);
      info.cput_lop = info.cput_lop+toc(time);
      if info.flag; return; end

      if ~isempty(d)	% the linear inequality system is feasible
        if nvp == p; info.ns = info.ns+1; end
        if verb2
          [info] = isf_print(V,bvec,perm,'lp',d,info,options,values);
          if info.flag; return; end
        end
        if p > nvp
          [info] = isf_rec(V,bvec,perm,d,info,options,values);
          if info.flag; return; end
        else
          if options.s
            if info.ns > info.nsmax
              info.s     = [info.s;zeros(info.nsmax,p)];
              info.nsmax = info.nsmax*2;
            end
            info.s(info.ns,perm) = bvec';
          end
        end
        info.nb_npl(nvp) = info.nb_npl(nvp)+1;
      else		% the linear inequality system is infeasible
        if verb2 >= 3
          fprintf(fout2,'\n%s | ',repmat(' ',1,floor(log10(info.ns))+1));
          [info] = isf_print_bvec(fout2,verb2,p,bvec,perm,info);
          fprintf(fout2,' | lp | infeasible system');
        end
        if options.sc; [info] = isf_memo_sc(p,nvp,bvec,perm,info,options); end
        if options.sv >= 2 	% add a new stem vector to info.stems
          [info] = isf_sv_add((2*bvec-1).*[lm;1],perm,info,options);	% update info.stems by adding a sign vector if appropriate
        end
      end

    end

  end

  return
