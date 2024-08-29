function [options,info] = isf_prelim(options,values)

%%
% [options,info] = isf_prelim(options,values)
%
% Do the preliminary verifications and decode options.

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

% Set default output values

  info.flag    = values.success;

% Check options fields

  possible_options = { ...
    'bestv', ...
    'dvnear0', ...
    'fout', ...
    'fout2', ...
    'prestring', ...
    'ranktol', ...
    'rc2018', ...
    's', ...
    'ssc', ...
    'sc', ...
    'sv', ...
    'svsprod', ...
    'verb', ...
    'verb2', ...
    'withd'};
  actual_options = fieldnames(options);
  for i=1:length(actual_options)
    if ~ismember(actual_options(i,:),possible_options)
      fprintf('\n\n### isf_prelim: field ''%s'' of structure ''options'' is not recognized\n\n',char(actual_options(i)));
      info.flag = values.fail_on_argument;
      return
    end
  end

% Set fout, verb, prestring

  if isfield(options,'verb')
    if isempty(options.verb) || (options.verb < 0) || (fix(options.verb) ~= options.verb)
      fprintf('\n\n### isf_prelim: options.verb must be a nonnegative integer\n\n');
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.verb = 1;			% default verbosity level
  end
  verb = options.verb;			% to reduce time access

  if isfield(options,'verb2')
    if isempty(options.verb2) || (options.verb2 < -1) || (fix(options.verb2) ~= options.verb2)
      fprintf('\n\n### isf_prelim: options.verb2 must be a nonnegative integer\n\n');
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.verb2 = 1;			% default verbosity level
  end

  if isfield(options,'fout')
    if isempty(options.fout) || (options.fout <= 0) || (fix(options.fout) ~= options.fout)
      if verb, fprintf('\n\n### isf_prelim: options.fout = %g not recongnized\n\n',options.fout); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.fout = 1;			% default output channel
  end
  fout = options.fout;			% to reduce time access

  if isfield(options,'fout2')
    if isempty(options.fout2) || (options.fout2 < 0) || (fix(options.fout2) ~= options.fout2)
      if verb, fprintf('\n\n### isf_prelim: options.fout2 = %g not recongnized\n\n',options.fout2); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.fout2 = 1;			% default output channel
  end

  if isfield(options,'prestring')
    if isempty(options.prestring) || ~ischar(options.prestring)
      fprintf('\n\n### isf_prelim: options.prestring must be a character array\n\n');
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.prestring = '';		% default verbosity level
  end
  prestring = options.prestring;	% to reduce time access

% Simplified algorithm of cerny-rada-2018 (rc2018)

  if isfield(options,'rc2018')
    if isempty(options.rc2018) || ~islogical(options.rc2018)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.rc2018 = %g must be a logical value\n\n',options.rc2018); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.rc2018 = false;	% default rc2018
  end

% Use modification B: if v'*d is near zero, determine 2 descendants

  if isfield(options,'dvnear0')
    if isempty(options.dvnear0) || ~islogical(options.dvnear0)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.dvnear0 = %g must be a logical value\n\n',options.dvnear0); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.dvnear0 = true;	% default dvnear0
  end

% Use modification D: construct in info.sv a list of stem vectors

  if isfield(options,'sv')
    if ~isnumeric(options.sv) || (fix(options.sv) ~= options.sv) || (options.sv < 0) || (options.sv > 3)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.sv (= %g) must be an integer in [0:3]\n\n',options.sv); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.sv = 3;		% default sv
  end

% Specify whether d must be computed at each node of the S-tree (this can be avoided only if options.sv == 3)

  if isfield(options,'withd')
    if isempty(options.withd) || ~islogical(options.withd)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.withd = %g must be a logical value\n\n',options.withd); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    if options.sv == 3
      options.withd = false;	% default withd
    else
      options.withd = true;	% default withd
    end
  end

% Use modification D: construct in info.svsprod a list of stem vectors

  if isfield(options,'svsprod')
    if ~isnumeric(options.svsprod) || (fix(options.svsprod) ~= options.svsprod) || (options.svsprod < 0) || (options.svsprod > 1)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.svsprod (= %g) must be an integer in [0:1]\n\n',options.svsprod); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.svsprod = 0;		% default svsprod
  end

% Use modification C: Take for the next v the best one in a sense in [1:3]; if 0, the next v index is the one in the 'perm' list

  if isfield(options,'bestv')
    if ~isnumeric(options.bestv) || (fix(options.bestv) ~= options.bestv) || ( (options.bestv ~= 0) && (options.bestv ~= 3) )
      if verb, fprintf(fout,'\n\n### isf_prelim: options.bestv (= %g) must be an integer in {0,3}\n\n',options.bestv); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    if options.sv == 3
      options.bestv = 0;	% default bestv
    else
      options.bestv = 3;	% default bestv
    end
  end

% List the feasible s's in info.s (s)

  if isfield(options,'s')
    if isempty(options.s) || ~islogical(options.s)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.s = %g must be a logical value\n\n',options.s); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.s = true;		% default s
  end

% List the infeasible s's: list nothing (0, default), short list (1), long list (2)

  if isfield(options,'sc')
    if isempty(options.sc) || ~islogical(options.sc)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.sc = %g must be a logical value\n\n',options.sc); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.sc = false;		% default sc
  end

% Specify the need of computing S and/or Sc

  if isfield(options,'ssc')
    if ~isnumeric(options.ssc) || (fix(options.ssc) ~= options.ssc) || (options.ssc < 1) || (options.ssc > 3)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.ssc (= %g) must be an integer in [1:3]\n\n',options.ssc); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.ssc = 1;		% default: compute S only
  end

% Specify the relative tolerance for computing the rank by the QR factorization

  if isfield(options,'ranktol')
    if ~isnumeric(options.ranktol) || (options.ranktol <= 0)
      if verb, fprintf(fout,'\n\n### isf_prelim: options.ranktol (= %g) must be a positive number\n\n',options.ranktol); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.ranktol = 1.e2*eps;	% default rank relative tolerance
  end

%-----------------------------------------------------------------------
% Check compatibility of the options
%-----------------------------------------------------------------------

  if options.sv == 3

    if (options.bestv > 0) && ~options.withd
      if verb, fprintf(fout,'\n\n### isf_prelim: options.bestv must be set to 0 when options.sv == 3 and options.withd is false\n\n'); end
      info.flag = values.fail_on_argument;
      return
    end

  else

    if options.ssc > 1
      if verb, fprintf(fout,'\n\n### isf_prelim: options.ssc (= %i) but the code for computing Sc is not available when options.sv ~= 3\n\n',options.ssc); end
      info.flag = values.fail_on_argument;
      return
    end

    if ~options.withd
      if verb, fprintf(fout,'\n\n### isf_prelim: options.withd (= %s) must be true when options.sv ~= 3\n\n',mat2str(options.withd)); end
      info.flag = values.fail_on_argument;
      return
    end

  end

  if (options.verb2 == 4) && ~options.withd
    if verb, fprintf(fout,'\n\n### isf_prelim: options.verb2 = 4 is not compatible with options.withd = false\n\n'); end
    info.flag = values.fail_on_argument;
    return
  end

  if (options.svsprod > 0) && (options.sv == 2)
    if verb, fprintf(fout,'\n\n### isf_prelim: options.svsprod = %i cannot be used when options.sv = 2\n\n',options.svsprod); end
    info.flag = values.fail_on_argument;
    return
  end

return
