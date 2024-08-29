function [info] = bf(V,options)

%%
% [info] = bf(V,options)
%
% Brute force algorithm, which list the sign vectors s such that
% s.*(V'*d) > 0 is feasible for d, using an optimization solver for min
% {t
% For more information, see section 5.2.1 in the paper 'On the
% B-differential of the componentwise minimum of two affine vector
% functions -- The full report' by J.-P. Dussault, J.Ch. Gilbert and B.
% Plaqevent-Jourdain [hal-03872711].
%
% On return, 'info' is a structure giving information on the run.
% . info.ns = half the number of feasible sign vectors;
% . info.s = half of the feasible sign vectors if options.s

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

% Set diagnosis values

  values.success                   = int8( 0);	% the required job has been realized

% Set default output variables

  info.flag = values.success;

% Get options

  fout      = options.fout;
  fout2     = options.fout2;
  verb      = options.verb;
  verb2     = options.verb2;

% Get dimensions

  [n,p] = size(V);

% Start printing

  if verb >= 2
    options.dline  = '------------------------------------------------------------------------';
    options.eline  = '========================================================================';
    options.tline  = '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~';
    fprintf(fout,'\n%s',options.eline);
  end
  if verb >= 1
    current_date = fix(datevec(now));
    fprintf(fout,'\nBF (Version 1.1, 2023-09-01)');
    fprintf(fout,'\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^');
    fprintf(fout,'\nCurrent date (time): %i-%02i-%02i (%i:%i)\n',current_date(1:5));
    fprintf(fout,'\nVector dimension (n)  = %i',n);
    fprintf(fout,'\nNumber of vectors (p) = %i',p);
    if verb2 == 1
      fprintf(fout,'\n\nBinary\nvector\n');
    elseif verb2 >= 2
      fprintf(fout,'\n\n# | Binary vector | LO optimal value | Direction\n\n');
    else
      fprintf(fout,'\n\n');
    end
  end

% Run the brute force method

  nsmax  = 2*p;		% max number of sign vector (will be increased below)
  info.s = zeros(nsmax,p);

  % Loop on the 2^p nodes. bvec = (s+1)/2 = binary representation of the sign vector s = (2*bvec)-1

  bvec = ones(p,1);	% initial binary vector made of all ones (any sign vector would be fine, it is a choice that allows us to list them all)
  ns = 0;		% start without selected signed vector

  % Main loop on the half of the binary vectors bvec making (2*bvec-1).*(V'*d) > 0 feasible for d. One always has bvec(1) == 1,
  % so that only half of the possible binary vectors are considered.

  while 1

    % Check wheter the sign vector yields a feasible system

    VS  = (2*bvec-1)'.*V;
    [d,val] = bf_optim(VS);

    if val <0
      ns = ns+1;
      if ns > nsmax
        info.s = [info.s;zeros(nsmax,p)];
        nsmax  = nsmax*2;
      end
      info.s(ns,:) = bvec';
      if verb2 == 1
        fprintf(fout2,'\n');
        fprintf(fout2,'%0i',bvec);
      end
      if verb2 == 2
        bf_print(bvec,d,val,ns,options);
      end
    end

    % Printing

    if verb2 >= 3
      bf_print(bvec,d,val,ns,options);
    end

    % Next binary vector

    b = bf_bin_minus(bvec(2:p));
    if isempty(b); break; end
    bvec(2:p) = b;

  end

% Update info

  if ns > 0
    info.s = info.s(1:ns,:);	% purge the tail of info.s
  end
  info.ns = ns;

% Last printings

  if verb >= 1
    fprintf(fout,'\n\nNumber of sign vectors = %i',ns*2);
    elapsed_time = toc(time0);
    fprintf(fout,'\n\nElapsed time: %g sec',elapsed_time);
    if ns ~= 0
      fprintf(fout,' (per sign vector: %g sec)',elapsed_time/(ns*2));
    end
    fprintf(fout,'\n%s',options.eline);
  end

return
