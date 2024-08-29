function [info] = isf_somesv(V,perm,info,options)

%%
% [info] = isf_somesv(V,perm,info,options)
%
% Compute the stem vectors of V that can be deduced from the r linear
% independent vectors given by the QR factorization and put them in the
% matrix info.stems of size [info.nb_stems,p].

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

% Get options

  fout      = options.fout;
  verb2     = options.verb2;

% Dimensions

  p = size(V,2);
  r = info.r;

% Precision parameter

  eps5 = 1.e5*eps;	% small positive number to detect nonzero elements in the computed null space vecctor

% Initialization

  info.nb_stems            = p-r;	% counter of stem vectors
  info.nb_duplicated_stems = 0;		% nb of computed stem vectors that are already in the matrix info.stems
  info.max_stems           = p-r;	% max nb of stem vectors
  info.stems               = zeros(info.max_stems,p);	% stem vector table in the form of a matrix

%-------------------------------------------------------------------------------------------------------------------------------
% Compute the p-r stem vectors deduced from the r linear independent vectors given by the QR factorization
%-------------------------------------------------------------------------------------------------------------------------------

  time = tic;

  i = 0;	% index of the stem vector

  for ivp = perm(r+1:p)		% index for a possible (r+1)th vector (perm is indeed horizontal)

    % construct a potentially new stem vector in 'stem'

    cols          = [perm(1:r),ivp];		% column indices of the considered vectors (horizontal)
    Z             = null(V(:,cols));		% it is a vector since nullity(V(:,cols)) == 1 (vertical)
    I             = find(abs(Z)>eps5);		% nonzero components of the null space vector
    stem          = zeros(1,p);
    stem(cols(I)) = sign(Z(I));			% one of the two stem vectors

    % memorize the stem vector

    i               = i+1;
    info.stems(i,:) = stem;

  end

%-------------------------------------------------------------------------------------------------------------------------------
% Purge info.stems of its duplicates (is this necesary in this very special case?)
%-------------------------------------------------------------------------------------------------------------------------------

  info.stems = unique(info.stems(1:info.nb_stems,:),'rows');

  nb_stems                 = size(info.stems,1);
  info.nb_duplicated_stems = info.nb_stems - nb_stems;
  info.nb_stems            = nb_stems;

%-------------------------------------------------------------------------------------------------------------------------------
% Compute the number of signs (nonzero components) for each stem vector, which may differ from r
%-------------------------------------------------------------------------------------------------------------------------------

  info.stem_sizes = sum(abs(info.stems),2);

%-------------------------------------------------------------------------------------------------------------------------------
% Printings
%-------------------------------------------------------------------------------------------------------------------------------

  if verb2 >= 2
    fprintf(fout,'\nNumber of stem vectors     = %i',info.nb_stems*2);
    fprintf(fout,'\nNumber of duplicated stems = %i',info.nb_duplicated_stems*2);
    fprintf(fout,'\nComputing time             = %g sec',toc(time));
    fprintf(fout,'\nStem vectors (without their symmetric):');
    for i=1:info.nb_stems
      stem = repmat(' ',1,p);
      stem(info.stems(i,:)==1) ='1';
      stem(info.stems(i,:)==-1)='0';
      fprintf('\n  | %s | (%i)',stem,i);
    end
    fprintf(fout,'\n%s',options.dline);
  end

%-------------------------------------------------------------------------------------------------------------------------------
% Compute the indices of the rows of info.stems with trailing zeros (trailing for the perm order). It is interesting to do this
% job once for all. This job is useless if options.bestv > 0 since the order of the vectors in perm can be modified at each
% node.
%-------------------------------------------------------------------------------------------------------------------------------

  if (options.bestv == 0) && (p > r+1)

    % info.stem_zero_indices{i} are the indices of the rows of info.stems with zeros in positions perm([p+1-i:p])

    info.stem_zero_indices = cell(p-r-1,1);

    % start with the determination of info.stem_zero_indices{1}

    info.stem_zero_indices{1} = find(info.stems(1:info.nb_stems,perm(p))==0);

    % determine info.stem_zero_indices{i} as a subset of info.stem_zero_indices{i-1}

    for i = 2:p-r-1
      I = info.stem_zero_indices{i-1};
      info.stem_zero_indices{i} = I(info.stems(I,perm(p+1-i))==0);
    end

  else

    info.stem_zero_indices = [];

  end

%-------------------------------------------------------------------------------------------------------------------------------

  return
