function [info] = isf_first(V,info,options,values)

%%
% [info] = isf_first(V,info,options,values)
%
% Prepare the recursif call to isf_rec.
%
% On return
% - info.ns: half the number of admissible sign vectors
%

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

%%%%%%%%%%%%%%%%%%%%%%
%% Table of Colors   %
%%%%%%%%%%%%%%%%%%%%%%
%% Bold   \033[1;--m %
%%%%%%%%%%%%%%%%%%%%%%
%% Black  \033[0;30m %
%% Red    \033[0;31m %
%% Green  \033[0;32m %
%% Yellow \033[0;33m %
%% Blue   \033[0;34m %
%% Purple \033[0;35m %
%% Cyan   \033[0;36m %
%% White  \033[0;37m %
%%%%%%%%%%%%%%%%%%%%%%

% Get options

  fout      = options.fout;
  verb      = options.verb;
  verb2     = options.verb2;

% Get dimensions

  [n,p] = size(V);

% Set default output values

  info.ns           = 0;		% nb of feasible sign vectors
  info.nsc          = 0;		% nb of infeasible sign vectors
  info.nb_losolve   = 0;		% nb of linear optimization problems solved
  info.nb_feaslop   = 0;		% nb of linear optimization problems yielding a feasible sign vector
  info.nb_infeaslop = 0;		% nb of linear optimization problems yielding an infeasible sign vector
  info.nb_svdetect  = 0;		% nb of infeasible s detected by a stem vector
  info.nb_new_stems = 0;		% nb of stem vectors discovered during the code running (when options.sv == 2)
  info.cput_lop     = 0;		% CPU time spent in solving LPs
  info.cput_sv      = 0;		% CPU time spent in computing stem vectors
  info.cput_cover   = 0;		% CPU time spent in checking that a sign vector covers a stem vector (out of the stem vector computation)

  info.flag         = values.success;

  info.counter.type2 = 0;	% number of type2 biconvex bipartitions
  info.counter.type3 = 0;	% number of type2, for which the separation failed, probably due to rounding errors
  info.counter.pert  = 0;	% number of handle perturbations
  info.counter.qp    = 0;	% number of QP solves

% Precision parameter

  eps2 = 1.e2*eps;	% small positive epsilon

% Preliminary printings

  if verb
    fprintf(fout,'\n%s',options.tline);
  end

%-------------------------------------------------------------------------------------------------------------------------------
% Normalize the vectors
%-------------------------------------------------------------------------------------------------------------------------------

  Vn2 = sum(V.*V,1); 	% row vector of the square norms of the V columns
  if any(Vn2==0)
    I = find(Vn2==0);
    if options.verb; fprintf(options.fout,'\n\n### isf_first: V(:,%i) is zero, nothing to do\n',I(1)); end
    info.flag = values.fail_on_zero_vector;
    return
  end
  V = V./sqrt(Vn2); 	% normalized vectors

%-------------------------------------------------------------------------------------------------------------------------------
% Detect the vectors that are 2 by 2 noncolinear and change V (the previous V is no longer necessary)
%-------------------------------------------------------------------------------------------------------------------------------

  p0 = p;
  [V,colsel,colout,info] = isf_noncolin(V,info,values);	% one can get rid of the original V
  if info.flag; return; end

  p = size(V,2);		% number of noncolinear vectors

  % Printings

  if verb >= 2
    if p ~= p0
      fprintf('\033[0;34m');
      fprintf(fout,'\nColinear vector(s):');
      eliminated_columns = setdiff(1:p0,colsel);
      fprintf(fout,'   %i//%i',[eliminated_columns;colout(eliminated_columns)]);
      fprintf('\033[0m');
      fprintf(fout,'\nVector dimension (n)  = %i',n);
      fprintf(fout,'\nNumber of vectors (p) = %i',p);
%     fprintf(fout,'\n%s',options.dline);
    end
  end

%-------------------------------------------------------------------------------------------------------------------------------
% Rank computation, using a QR factorization
%-------------------------------------------------------------------------------------------------------------------------------

  [Q,R,perm] = qr(V,0);	% V(:,perm) = V*P = Q*R and abs(diag(R)) is decreasing (QR should have "economy size", but it hasn't)

  Rnorm = norm(R,Inf);
  tol   = options.ranktol*Rnorm;
  r     = size(R,1);	% detect rank, up to options.ranktol
  while (r > 0) && (norm(R(r,:),1) < tol)
    r = r-1;
  end
  info.r = r;

  % Set default output values depending on p

  info.nb_npl   = zeros(p,1);	% nb of feasible sign vectors per level of the S-tree

  info.schlafli = isf_schlafli(p,info.r);
  nsmax         = info.schlafli/2;
  info.nsmax    = nsmax;		% length of info.s
  info.s        = zeros(nsmax,p);	% list of sign vectors

  if options.sc && (options.sv ~= 3)
    nscmax      = 2*p;
    info.nscmax = nscmax;		% length of info.sc
    info.sc     = zeros(nscmax,p);	% list of "infeasible" sign vectors, if this is required by the user
  end
  info.sv       = {};		% "infeasible" sign vectors used by the algorithm to prevent exploring parts of the S-tree

  % Set the number of nodes in the S-tree for the levels [1:r]

  info.nb_npl(1) = 1;
  for l = 2:r
    info.nb_npl(l) = info.nb_npl(l-1)*2;
  end

  if verb2 >= 2
    fprintf(fout,'\nRank (r)           = %i',r);
    fprintf(fout,'\nnorm(R,Inf)        = %8.2e',Rnorm);
    fprintf(fout,'\noptions.ranktol    = %8.2e',options.ranktol);
    fprintf(fout,'\nabsolute tolerance = %8.2e',tol);
    if r < p
      fprintf(fout,'\n%s\nSelected linearly independent vectors %i',options.dline,perm(1));
      fprintf(fout,', %i',perm(2:r));
      fprintf(fout,'\nOther linearly dependent vectors      %i',perm(r+1));
      for i = r+2:p; fprintf(fout,', %i',perm(i)); end
      fprintf(fout,'\n%s',options.dline);
    end
  end

%-------------------------------------------------------------------------------------------------------------------------------
% Special case when all the stem vectors must be compted
%-------------------------------------------------------------------------------------------------------------------------------

  if options.sv == 3

    time = tic;
    [info] = isf_sv_all(V,perm,info,options,values);	% info.stems does not take perm into account
    info.cput_sv = toc(time);

    if (options.ssc == 2) || (options.ssc == 3)		% the computation of Sc is required
      fprintf('\n\n#### isf_first: the computation of Sc has not yet been implemented\n');
      fprintf(fout,'\n%s',options.dline);
    end

    if options.ssc == 2		% the computation of S is not required
      return
    end

%-------------------------------------------------------------------------------------------------------------------------------
% Compute the 2(p-r) stem vectors that can be deduced from the r linear independent vectors given by the QR factorization and
% put them in info.stems
%-------------------------------------------------------------------------------------------------------------------------------

  elseif options.sv > 0

    if p == r

      info.stems      = [];
      info.nb_stems   = 0;
      info.cput_sv    = 0;
      info.cput_cover = 0;

    else

      time = tic;
      [info] = isf_somesv(V,perm,info,options);
      info.cput_sv = toc(time);

    end

  end

%-------------------------------------------------------------------------------------------------------------------------------
% Run the incremental method
%-------------------------------------------------------------------------------------------------------------------------------

  % The vector bvec = (s+1)/2 is the binary representation of the sign vector s = (2*bvec)-1. It is linked to the current node
  % and its length is the number of vectors considered at the current node. It applies to the column vectors of
  % V(:,perm(1:length(bvec))).

  bvec = ones(r,1);	% initial binary vector made of all ones (any sign vector would be fine, it is a choice that allows us to list them all)

  RmT = inv(R(1:r,1:r))';	% there are 2^(r-1) linear systems to solve with R', which is more than r=size(R,1) when r ≥ 3

  % Main loop on the half of the binary vectors bvec making (2*bvec-1).*(V'*d) > 0 feasible for d. One always has bvec(1) == 1,
  % so that only half of the possible binary vectors are considered.

  if ~options.withd; d = []; end

  while 1

    % Product of the stem-vectors times the current sign vector

    if options.svsprod == 1
      time = tic;
      s               = 2*bvec-1;					% sign vector of the current node of the S-tree (is generated inpendently of perm)
      info.svsprod    = info.stems(1:info.nb_stems,perm(1:r))*s;	% value transmitted to the next S-tree level
      info.cput_cover = info.cput_cover+toc(time);
    end

    % Associated direction

    if options.withd
      d = Q(:,1:r)*(RmT*(2*bvec-1));
      d = d/norm(d);
    end

    % Printing

    if verb2
      if (verb2 >= 3) || (r == p)
        [info] = isf_print(V,bvec,perm,'qr',d,info,options,values);
        if info.flag; return; end
      end
    end

    % recursive calls

    if r < p
      if options.withd
        [info] = isf_rec(V,bvec,perm,d,info,options,values);
      else
        [info] = isf_rec_nod(V,bvec,perm,info,options,values);
      end
      if info.flag
        return
      end
    else
      info.ns = info.ns+1;
      info.s(info.ns,perm) = bvec';
    end

    % Next binary vector

    b = isf_bin_minus(bvec(2:r));
    if isempty(b); break; end
    bvec(2:r) = b;		% maintain bvec(1) = 1 (only half the S-tree needs to be generated, by symmetry)

  end

%-------------------------------------------------------------------------------------------------------------------------------
% In case there were colinear columns in the original V, take this into account to update info.s and info.sc
%-------------------------------------------------------------------------------------------------------------------------------

  if p ~= p0

    s      = info.s;			% sign vectors in binary representation of the "transformed V"
    info.s = zeros(info.ns,p0);		% aims at being the sign vectors in binary representation of the "original V"
    info.s(:,colsel) = s;
    cols = eliminated_columns(colout(eliminated_columns)>0);
    cols_to_copy = colout(cols);
    info.s(:,cols) = info.s(:,cols_to_copy);
    cols = eliminated_columns(colout(eliminated_columns)<0);
    cols_to_copy = -colout(cols);
    info.s(:,cols) = 1-info.s(:,cols_to_copy);

    if verb >= 2
      fprintf(fout,'\n%s',options.tline);
      fprintf('\033[0;34m');
      fprintf(fout,'\nBack to the original %i vector setting',p0);
      fprintf('\033[0m');
      print_table(info.s,options);
    end

  end

%-------------------------------------------------------------------------------------------------------------------------------
% Return
%-------------------------------------------------------------------------------------------------------------------------------

  if verb
    fprintf(fout,'\n%s',options.tline);
  end

  return

%%%%%%%%%%%%%%%%%%%%%
%% Nested function %%
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = print_table(blist,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  fout      = options.fout;
  prestring = options.prestring;

  fprintf(fout,'\n%sDirect and complementary sign vectors',prestring);

  for i = 1:size(blist,1)
    fprintf(fout,'\n%s  %2i:  ',prestring,i);
    fprintf(fout,'%i',blist(i,:));
    fprintf(fout,' | ');
    fprintf(fout,'%i',1-blist(i,:));
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% End of nested function %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
