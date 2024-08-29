function [info] = isf(V,options)

%%
% [info] = isf(V,options)
%
% For p real nonzero vectors of dimension n, stored in the columns of
% the nxp matrix V, denote by S the set of sign vetors s in {-1,+1}^p
% such that s.*(V'*d)>0 is feasible for d in R^n and denote by S^c =
% {-1,+1}^p\S the complementary set of S. ISF computes half the sign
% vectors in S (and in S^c, on request). The other half can be obtained
% by symmetry (see below). The selected sign vectors are the leaves of a
% tree of sign vectors of increasing size, called the S-tree below.
%
% The 'binary representation' of a sign vector s in {-1,+1}^p is
% (s+1)/2, which is indeed only formed of elements in {0,1}. This binary
% representation is useful, in particular, for printing s in a compact
% manner.
%
% ISF stands for Incremental Sign Feasibility (not 'Impôt Sur la
% Fortune'). The term `Incremental' refers to the fact that the
% algorithm constructs the S-tree of the sign vectors incrementally. The
% term 'Sign' refers to the sign vectors. The term 'Feasibility' refers
% to the fact that the feasibility of s.*(V'*d) > 0 in d is used to
% select the appropriate sign vectors. Finally, the juxtaposition of
% terms in 'Incremental Sign Feasibility' is meaningless.
%
% The code is described in the paper 'On the B-differential of the
% componentwise minimum of two affine vector functions' by J.-P.
% Dussault, J.Ch. Gilbert and B. Plaqevent-Jourdain. This paper
% submitted to Mathematical Programming Computation is available on
% http://hal.archives-ouvertes.fr/hal-04048393.
%
% Various options allow the user to tune the code. All options have
% default values. Below, LOP stands for 'linear optimization problem'.
% . options.bestv (integer in {0,3}; default 0 if options.sv == 3;
%   default 3 if options.sv < 3; ignored if options.rc2018 == true or
%   options.withd == false): can be used to modify the order in which
%   the vectors (columns of V) are considered for constructing the
%   S-tree; be default, this order is that imposed by the QR
%   factorization of V; only the value 0 and 3 are actually documented
%   and evaluated in the paper:
%   == 0: no modification of the order, except that imposed by the QR
%         factorization;
%   == 3: a reodering of the vectors is made so that the number of nodes
%         of the S-tree decreases (and therefore the number of LOPs to
%         solve);
% . options.dvnear0 (logical; default value is 'true'; ignored if
%   options.rc2018 == true or options.withd == false): if true, the
%   S-tree algorithm constructs two descendants when v'*d is in a
%   specific computed interval surrounding zero, without having to solve
%   a LOP or to use stem vectors if any (v is the new considered
%   vector);
% . options.fout: first output channel containing all but the printing
%   made during the generation of the S-tree; set 1 (default) for the
%   standard output (screen); for a specific file, use "options.fout =
%   fopen(...)" before calling ISF;
% . options.fout2: second output channel containing the printing made
%   during the generation of the S-tree; set 1 (default) for the
%   standard output (screen); for a specific file, use "options.fout2 =
%   fopen(...)" before calling ISF;
% . options.ranktol (real number, default '1.e2*eps'): gives the
%   relative tolerance for computing the rank of V by its QR
%   factorization; a row of R is considered to be zero if its l1-norm is
%   less than options.ranktol*norm(R,Inf);
% . options.rc2018 (logical, default 'false'): if true, the
%   simulated Rada and Černý algorithm is run, ignoring most other
%   options;
% . options.s (logical; default 'true'): if 'true', half of the sign
%   vectors s in S are stored in info.s; the other half can be obtained
%   by symmetry: ones(size(info.s))-info.s;
% . options.sc (logical; default value is 'false'): if 'true' and
%   options.sv < 3, the infeasible s's are stored in info.sc; the other
%   half can be obtained by symmetry: ones(size(info.sc))-info.sc;
% . options.sv (integer in [0:3]; default 3; ignored if options.rc2018
%   == true): determine whether and how many stem vectors must be used
%   (it is the option Di in the paper);
%   == 0: do not use stem vectors;
%   == 1: use the p-r stem vectors that can be computed thanks to the QR
%         factorization of V (r is the rank of V); this is not many, but
%         it is inexpensive to compute;
%   == 2: in addition to the p-r stem vectors obtained with
%         options.sv=1, use the dual solution of each LO problem having
%         no solution to add a stem vector to the list; this improves
%         significantly the speed of the algorithm, but requires to
%         solve the LOPs with the dual simplex method;
%   == 3: compute all the stem vectors at the begining of the run
%         (time consuming operation when p is large);
% . options.svsprod (integer in [0:1]; default 1 if (options.sv == 3) &&
%   ~options.withd, 0 otherwise; must be 0 if options.sv == 2) determine
%   how is detected the fact that a sign vector covers a stem vector;
%   == 0: the detection is done in the function isf_cover_stem by making
%         a product of the matrix of the stem vectors and the sign
%         vector;
%   == 1: partial dot products are generated at each node of the
%         S-tree;
% . options.verb (a nonnegative number, default 1): verbosity level of
%   the output channel options.fout:
%   == 0: works silently;
%   >= 1: error messages, initial setting, final status (default);
%   >= 2: more information;
% . options.verb2: verbosity level of the output on channel
%   options.fout2;
%   == 0: works silently;
%   >= 1: error messages and binary representation of the sign vectors;
%   >= 2: the feasible directions d;
%   >= 3: some information at the intermediate steps of the recursivity
%         process;
%   >= 4: verification that the directions d verify s.*(V'*d) > 0;
% . options.withd (logical, default 'true' if options.sv < 3, default
%   'false' if options.sv == 3): if 'true', a direction d is
%   associated with (and computed at) each node of the S-tree; this
%   requires more computation, in particular for solving LOPs.
%
% On return, 'info' is a structure giving information on the run.
% . info.flag (integer in [0:9]): diagonis of the run;
%   == 0: the required job has been realized;
%   == 1: an input argument is wrong;
%   == 2: nothing to do since V has no column (p=0),
%   == 3: one of the vectors V(:,i) vanishes (then no solution),
%   == 4: at some node of the S-tree, a computed d does not satisfy
%         z.*(W'*d)>0, while it should (here z is a sign vector that may
%         have less that p components and W may be formed of part of the
%         columns of V); this is usually due to rounding error, although
%         an implementation error cannot be excluded;
%   == 5: one vector is opposite to another one, in which case the
%         problem has no solution;
%   == 6: the linear optimization solver failed;
%   == 9: a technical problem is encountered, which requires
%         improvements of the code;
% . info.ns = half the number of feasible sign vectors in S;
% . info.nsc = half the number of infeasible sign vectors in S^c;
% . info.s = matrix of size [info.ns,p] containing half of the feasible
%   sign vectors if options.s == true; the other half can be obtained by
%   symmetry: ones(size(info.s))-info.s;
% . info.sc = (if options.sc == true and options.sv < 3) matrix of size
%   [info.nsc,p] containing half of the infeasible sign vectors in S^c;
%   the other half can be obtained by symmetry:
%   ones(size(info.sc))-info.sc.

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

% Start timer

  time0 = tic;

% Set default output values

  info = [];

% Set diagnosis values

  values.success                   = int8( 0);	% the required job has been realized
  values.fail_on_argument          = int8( 1);	% an input argument is wrong
  values.nothing_to_do             = int8( 2);	% nothing to do since V has no column
  values.fail_on_zero_vector       = int8( 3);	% one of the vectors V(:,i) vanishes (then no solution)
  values.fail_on_pointed_cone      = int8( 4);	% a pointed cone and its handle are not compatible
  values.fail_on_opposite_vector   = int8( 5);	% one vector is opposite to another one
  values.fail_on_lp_solve          = int8( 6);	% the LP solver fails
  values.fail_on_special_lp        = int8( 7);	% unexpected case with linprog, requires more attention
  values.fail_on_qp_solve          = int8( 8);	% the QP solver fails
  values.fail_on_sd_not_compatible = int8( 9);	% a s(i)*v(i)'*d <= for some i
  values.fail_on_technicality      = int8(10);	% a technical problem is encountered, which requires improvements

% Get options

  if (nargin < 2) || isempty(options)
    options = struct();
  end
  if ~isstruct(options)
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','Argument ''options'' is expected to be a structure\n');
  end

  [options,info] = isf_prelim(options,values);
  if info.flag; return; end

  % to reduce time access, set:

  fout      = options.fout;
  verb      = options.verb;
  prestring = options.prestring;

% Get dimensions

  [n,p] = size(V);

  if p == 0
    info.ns = 0;
    info.flag = values.nothing_to_do;
    if verb, fprintf(fout,'\n\n### isf: nothing to do since V has no column\n\n'); end
    return
  end

% Start printing

  if verb
    options.dline  = '------------------------------------------------------------------------';
    options.eline  = '========================================================================';
    options.tline  = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
    fprintf(fout,'\n%s%s',prestring,options.eline);

    current_date = fix(datevec(now));
    fprintf(fout,'\n%sISF (Version 1.1, 2023-09-01)',prestring);
    fprintf(fout,'\n%s^^^^^^^^^^^^^^^^^^^^^^^^^^^^^',prestring);
    fprintf(fout,'\n%sCurrent date (time): %i-%02i-%02i (%i:%i)\n',prestring,current_date(1:5));
    fprintf(fout,'\n%sVector dimension (n)  = %i',prestring,n);
    fprintf(fout,'\n%sNumber of vectors (p) = %i',prestring,p);
    fprintf(fout,'\nAlgorithm features');
    fprintf(fout,'\n. Incremental-recursive method');
    if options.rc2018
      fprintf(fout,'\n. Rada and Černý (2018) simulated version');
    end
    if options.withd
      fprintf(fout,'\n. A direction d is computed at each node of the S-tree');
    else
      fprintf(fout,'\n. No direction d is computed');
    end
    if (options.sv ~= 3) || ( (options.sv == 3) && options.withd )
      if options.dvnear0
        fprintf(fout,'\n. Use v''*d near zero, to determine 2 descendants (option B)');
      else
        fprintf(fout,'\n. Use v''*d=0, to determine 2 descendants');
      end
    end
    if options.bestv > 0
      fprintf(fout,'\n. The next v is the best one in the sense C%0i',options.bestv);
    else
      fprintf(fout,'\n. The next v is the one in the list established by the QR fctorization with permutations (option C0)');
    end
    if options.sv == 1
      fprintf(fout,'\n. Construct a list of stem vectors after the QR factorization');
    elseif options.sv == 2
      fprintf(fout,'\n. Construct a list of stem vectors after the QR factorization and with the LO solved for infeasible s');
    elseif options.sv == 3
      fprintf(fout,'\n. Construct the list of all the stem vectors');
      if options.svsprod > 0
        fprintf(fout,'\n. Matrix-vector product for determining whether a sign vector covers a stem vector is scattered in the S-tree');
      end
    end
    if options.ssc == 1
      fprintf(fout,'\n. Compute S');
    elseif options.ssc == 2
      fprintf(fout,'\n. Compute Sc');
    else
      fprintf(fout,'\n. Compute S and Sc');
    end
    fprintf(fout,'\nOutput');
    fprintf(fout,'\n. channel 1 (fout)          = %0i',options.fout);
    fprintf(fout,'\n. verbosity level 1 (verb)  = %0i',options.verb);
    if options.verb2 >= 0
      fprintf(fout,'\n. channel 2 (fout2)         = %0i',options.fout2);
    end
    fprintf(fout,'\n. verbosity level 2 (verb2) = %0i',options.verb2);
  end

% Run the incremental method

  if options.rc2018
    [info] = isf_first_rc2018(V,info,options,values);
  else
    [info] = isf_first(V,info,options,values);
  end

  info.cput_total = toc(time0);

% Update info

  if (options.s) && (info.ns > 0)
    info.s = info.s(1:info.ns,:);	% purge the tail of info.s
  end

  if options.sc
    info.sc = info.sc(1:info.nsc,:);	% purge the tail of info.sc
  else
    info.nsc = 2^(p-1)-info.ns;		% = |S^c|/2
  end

% Last printings

  if verb
    if info.flag
      if info.flag == values.fail_on_lp_solve
        fprintf('\n### Linprog failed to solve a problem');
      elseif info.flag == values.fail_on_special_lp
        fprintf('\n### Unexpected case with Linprog, which requires more attention');
      end
    else
      fprintf(fout,'\nStatistics');
      fprintf(fout,'\n. number of sign vectors           = %i (≤ %0i = Schläfli''s upper bound)',info.ns*2,info.schlafli);
      if options.sv == 3
        if options.withd; fprintf(fout,'\n. number of LO problems solved     = %i (feas %0i, infeas %0i)',info.nb_losolve,info.nb_feaslop,info.nb_infeaslop); end
      else
        fprintf(fout,'\n. number of LO problems solved     = %i (feas %0i, infeas %0i)',info.nb_losolve,info.nb_feaslop,info.nb_infeaslop);
      end
      if options.sv > 0 && ~options.rc2018
        fprintf(fout,'\n. number of initial stem vectors   = %i',info.nb_stems-info.nb_new_stems);
        if options.sv == 2
          fprintf(fout,'\n. number of running stem vectors   = %i',info.nb_new_stems);
        end
        fprintf(fout,'\n. number of stem vector detections = %i',info.nb_svdetect);
      end
      fprintf(fout,'\nS-tree (half of it)\n. Nodes per level:');
      fprintf(fout,' %0i',info.nb_npl);
      fprintf(fout,'\n. Total number of nodes: %0i',sum(info.nb_npl));
      fprintf(fout,'\nCPU time');
      if options.sv > 0 && ~options.rc2018
        fprintf(fout,'\n. initial stem vectors  %9.2e sec (%2.1f %%)',info.cput_sv,info.cput_sv/info.cput_total*100);
      end
      if (options.sv < 3) || options.withd
        fprintf(fout,'\n. LO solves             %9.2e sec (%2.1f %%)',info.cput_lop,info.cput_lop/info.cput_total*100);
      end
      if options.sv > 0 && ~options.rc2018
        fprintf(fout,'\n. cover detection       %9.2e sec (%2.1f %%)',info.cput_cover,info.cput_cover/info.cput_total*100);
      end
      fprintf(fout,'\n. total                 %9.2e sec',info.cput_total);
      if info.ns ~= 0
        fprintf(fout,'\n. per sign vector       %9.2e sec',info.cput_total/(info.ns*2));
      end
    end
  end
  if verb >= 2
    fprintf(fout,'\n%s%s\n',prestring,options.eline);
  end

  return
