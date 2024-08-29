function [info] = isf_first_rc2018(V,info,options,values)

%%
% [info] = isf_first_rc2018(V,info,options,values)
%
% Prepare the recursif call to isf_rec_rc2018.
%
% On return
% - info.ns: half the number of admissible sign vectors

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

  nsmax = 2*p;

  info.ns           = 0;		% nb of feasible sign vectors
  info.nsmax        = nsmax;		% length of info.s
  info.s            = zeros(nsmax,p);	% list of sign vectors
  info.nb_npl       = zeros(p,1);	% nb of feasible sign vectors per level
  info.nb_losolve   = 0;		% nb of linear optimization problems solved
  info.nb_feaslop   = 0;		% nb of feasible linear optimization problems
  info.nb_infeaslop = 0;		% nb of infeasible linear optimization problems
  info.cput_lop     = 0;		% CPU time spent in solving LPs

  info.flag         = values.success;

  info.counter.type2 = 0;	% number of type2 biconvex bipartitions
  info.counter.type3 = 0;	% number of type2, for which the separation failed, probably due to rounding errors
  info.counter.pert  = 0;	% number of handle perturbations
  info.counter.qp    = 0;	% number of QP solves

% Normalize the vectors and detect those that are 2 by 2 noncolinear

  p0 = p;
  [V,colsel,colout,info] = isf_noncolin(V,info,values);	% one can get rid of the original V
  if info.flag; return; end

  p = size(V,2);
  if verb >= 2
    fprintf(fout,'\n%s',options.tline);
    if p ~= p0
      fprintf('\033[0;34m');
      fprintf(fout,'\nColinear vector(s):');
      eliminated_columns = setdiff(1:p0,colsel);
      fprintf(fout,'   %i//%i',[eliminated_columns;colout(eliminated_columns)]);
      fprintf('\033[0m');
      fprintf(fout,'\n%s',options.dline);
    end
  end

% Printings

  if verb >= 2
    fprintf(fout,'\nVector dimension (n)  = %i',n);
    fprintf(fout,'\nNumber of vectors (p) = %i',p);
%   fprintf(fout,'\n%s',options.dline);
  end

  % Set the number of nodes in the S-tree for the levels [1:r]

  info.nb_npl(1) = 1;

% Run the recursive method

  if verb >= 2
    fprintf(fout,'\n%s',options.dline);
  end

  % Call to the recursive process, which starts with bvec == 1,
  % so that only half of the possible binary vectors are considered.

  bvec = 1;	% initial binary vector made of all ones (any sign vector would be fine, it is a choice that allows us to list them all)

  % No permutation of the columns of V

  perm = 1:p;

  % Direction associated with s(1)=+1

  d = V(:,1);
  d = d/norm(d);

  % Printing

  if verb2 >= 3
    [info] = isf_print(V,bvec,perm,'  ',d,info,options,values);
    if info.flag; return; end
  end

  % Recursive call

  if p > 1
    [info] = isf_rec_rc2018(V,bvec,perm,d,info,options,values);
    if info.flag; return; end
  end

  if verb >= 2
    fprintf(fout,'\n%s',options.tline);
  end

return
