function [options,info] = bdiffmin_prelim(options,values)

%%
% [info] = bdiffmin_prelim(options,values)
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
% <http://doc.trolltech.com/3.0/license.html>.
%
%-----------------------------------------------------------------------

% Set default output values

  info.flag    = values.success;

% Check options fields

  possible_options = { ...
    'bdiffc', ...
    'eqtol', ...
    'fout', ...
    'verb'};
  actual_options = fieldnames(options);
  for i=1:length(actual_options)
    if ~ismember(actual_options(i,:),possible_options)
      fprintf('\n\n### bdiffmin_prelim: field ''%s'' of structure ''options'' is not recognized\n\n',char(actual_options(i)));
      info.flag = values.fail_on_argument;
      return
    end
  end

% Set fout, verb

  if isfield(options,'verb')
    if isempty(options.verb) || (options.verb < 0) || (fix(options.verb) ~= options.verb)
      fprintf('\n\n### bdiffmin_prelim: options.verb must be a nonnegative integer\n\n');
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.verb = 1;			% default verbosity level
  end
  verb = options.verb;			% to reduce time access

  if isfield(options,'fout')
    if isempty(options.fout) || (options.fout <= 0) || (fix(options.fout) ~= options.fout)
      if verb, fprintf('\n\n### bdiffmin_prelim: options.fout = %g not recongnized\n\n',options.fout); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.fout = 1;			% default output channel
  end
  fout = options.fout;			% to reduce time access

% Equality tolorence

  if isfield(options,'eqtol')
    if isempty(options.eqtol) || (options.eqtol <= 0)
      if verb, fprintf('\n\n### bdiffmin_prelim: options.eqtol must be a (small) positive number\n\n'); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.eqtol = 1.e-8;		% default equality tolerance
  end

% List of non Jacobians

  if isfield(options,'bdiffc')
    if isempty(options.bdiffc) || ~islogical(options.bdiffc)
      if verb, fprintf(fout,'\n\n### bdiffmin_prelim: options.bdiffc = %g must be a logical value\n\n',options.bdiffc); end
      info.flag = values.fail_on_argument;
      return
    end
  else
    options.bdiffc = false;		% default bdiffc
  end

return
