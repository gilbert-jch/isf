function [info] = bdiffmin(A,a,B,b,x,options)

%%
% [info] = bdiffmin(A,a,B,b,x,options)
%
% Bdiffmin list in the cell array in info.bdiff the elements of the
% B-differential of the function
%
%    H: x -> H(x) = min(A*x+a,B*x+b),
%
% where A and B are mxn matrices and a and b are mx1 vectors. Hence
% info.bdiff{i} is the i-th Jacobian of the B-differential of H at the
% given x, which is an mxn matrix. There are length(c) such Jacobians.
%
% The list of the Jacobians in the B-differential is obtained thanks to
% the function 'isf'.
%
% Various options allow the user to tune the code:
% . options.bdiffc (logical, default 'false'): if 'true' the Jacobians
%   that are in the Cartesian proiduct of the B-differential of the
%   components but not in the B-differential are also listed in
%   info.bdiffc;
% . options.eqtol: small value used to detect equality between the
%   components of A*x+a and B*x+b;
% . options.fout: output channel containing all the printings made
%   during the run of bdiffmin (nothing is printed from ISF); set 1
%   (default) for the standard output (screen); for a specific file, use
%   "options.fout = fopen(...)" before calling BDIFFMIN;
% . options.verb (a nonnegative number, default 1): verbosity level of
%   the output channel options.fout:
%   == 0: works silently;
%   >= 1: error messages, initial setting, final status (default).
%
% On return, info is a structure giving information on the run
% - info.bdiff = cell array. Its i-th element is the i-th Jacobian of
%   the B-differential of H. This is an nxn matrix.
% - info.bdiffc = cell array (if options.bdiffc == 'true'). Its i-th
%   element is the i-th matrix that is in the Cartesian product of the
%   B-differential of the components of H but not in the B-differential
%   of H;
% . info.flag (integer in [0:11]): diagonis of the run;
%   ==  0: the required job has been realized;
%   == 11: an input argument is wrong;
%   the other flags are inherited from isf.

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
% Bdiffmin is distributed under the terms of the Q Public License
% version 1.0.
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

% Set diagnosis values

  values.success                   = int8( 0);	% the required job has been realized
  values.fail_on_argument          = int8(11);	% an input argument is wrong

% Set default output values

  info       = [];
  info.flag  = values.success;

% Check input arguments

  if nargin < 5
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','There must be at least 5 input arguments\n');
  end
  [m,n] = size(A);
  if any([m,n]~=size(B))
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','The matrices ''A'' and ''B'' must have the same sizes\n');
  end
  if any([m,1]~=size(a))
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','The sizes of ''A'' and ''a'' are not compatible\n');
  end
  if any([m,1]~=size(b))
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','The sizes of ''B'' and ''b'' are not compatible\n');
  end
  if any([n,1]~=size(x))
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','The sizes of ''A'' and ''x'' are not compatible\n');
  end

% Decode options

  if (nargin < 6) || isempty(options)
    options = struct();
  end
  if ~isstruct(options)
    info.flag = values.fail_on_argument;
    error('isf:OptionsNotStruct','Argument ''options'' is expected to be a structure\n');	% followed by return
  end

  [options,info] = bdiffmin_prelim(options,values);
  if info.flag; return; end

  % to reduce time access

  eqtol = options.eqtol;
  fout  = options.fout;
  verb  = options.verb;

% Analyze data

  F = A*x+a;
  G = B*x+b;

  indE  = find(abs(G-F) <= eqtol);
  Idiff = setdiff(1:m,indE);
  indF  = Idiff(F(Idiff)<G(Idiff));
  indG  = setdiff(Idiff,indF);

  indV = zeros(m,1);	% indices of nondifferentiability
  j = 0;
  for i = indE(:)'
    if norm(B(i,:)-A(i,:),Inf) > eqtol
      j = j+1;
      indV(j) = i;
    end
  end
  indV = indV(1:j);
  indD = setdiff(1:m,indV);	% indices of differentiability
  V = (B(indV,:)-A(indV,:))';

% Start printing

  if verb

    options.eline  = '============================================================================';
    fprintf(fout,'\n%s',options.eline);

    current_date = fix(datevec(now));
    fprintf(fout,'\nBdiffmin (Version 1.1, 2023-09-01)');
    fprintf(fout,'\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^');
    fprintf(fout,'\nCurrent date (time): %i-%02i-%02i (%i:%i)\n',current_date(1:5));
    fprintf(fout,'\nNumber of functions (m) = %i',m);
    fprintf(fout,'\nNumber of variables (n) = %i',n);
    fprintf(fout,'\nTolerances');
    fprintf(fout,'\n. equality (eqtol) = %7.1e',eqtol);
    fprintf(fout,'\nComponent analysis');
    if ~isempty(indF)
      fprintf(fout,'\n. (A*x+a)(i) < (B*x+b)(i) for i =');
      fprintf(fout,' %i',indF);
    end
    if ~isempty(indG)
      fprintf(fout,'\n. (A*x+a)(i) > (B*x+b)(i) for i =');
      fprintf(fout,' %i',indG);
    end
%   if ~isempty(indE)
%     fprintf(fout,'\n. (A*x+a)(i) = (B*x+b)(i) for i =');
%     fprintf(fout,' %i',indE);
%   end
%   if ~isempty(indD)
%     fprintf(fout,'\n. min(A*x+a,B*x+b) is differentiable for i =');
%     fprintf(fout,' %i',indD);
%   end
    if ~isempty(setdiff(indE,indV))
      fprintf(fout,'\n. (A*x+a)(i) = (B*x+b)(i) and A(i,:) == B(i,:) for i =');
      fprintf(fout,' %i',setdiff(indE,indV));
    end
    if ~isempty(indV)
      fprintf(fout,'\n. (A*x+a)(i) = (B*x+b)(i) and A(i,:) ~= B(i,:) for i =');
      fprintf(fout,' %i',indV);
    end
    fprintf('\n\nRunning isf\n');
  end

% Run the incremental method

  options_isf.verb  = 0;
  options_isf.verb2 = 0;
  if options.bdiffc
    options_isf.bestv = 3;
    options_isf.sv    = 2;
    options_isf.sc    = true;
  end

  [info] = isf(V,options_isf);
  if info.flag; return; end

% Construct the B-differential

  % Common part of the Jacobians

  J         = zeros(m,n);
  J(indF,:) = A(indF,:);
  J(indD,:) = A(indD,:);
  J(indG,:) = B(indG,:);

  % Variable part of the Jacobians

  info.bdiff = {};
  blist = sortrows([info.s;ones(size(info.s))-info.s]);
  nj = 0;
  for i = 1:2*info.ns
    s = blist(i,:)';
    nj = nj+1;
    J(indV,:) = (1-s).*A(indV,:)+s.*B(indV,:);
    info.bdiff{nj} = J;
  end

% Construct the complement of the B-differential

  if options.bdiffc

    info.bdiffc = {};
    bclist = sortrows([info.sc;ones(size(info.sc))-info.sc]);
    nj = 0;
    for i = 1:2*info.nsc
      sc = bclist(i,:)';
      nj = nj+1;
      J(indV,:) = (1-sc).*A(indV,:)+sc.*B(indV,:);
      info.bdiffc{nj} = J;
    end

  end

% Final printings

  if verb
    fprintf(fout,'\nNumber of Jacobians: %i',info.ns*2);
    fprintf(fout,'\nNumber of non Jacobians: %i',info.nsc*2);
    fprintf(fout,'\n\nElapsed time: %g sec',toc(time0));
    fprintf(fout,'\n%s\n',options.eline);
  end

  return
