function [info] = isf_sv_all(V,perm,info,options,values)

%%
% [info] = isf_sv_all(V,perm,info,options,values)
%
% Compute all the stem vectors of V and put them in the matrix
% info.stems of size [info.nb_stems,p].

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
  if p < 3
    fprintf(fout,'\n\n#isf_sv_all: this function has not been designed for p < 3 (and p = %i)\n\n',p);
    info.flag = values.fail_on_technicality;
    return
  end
  r = info.r;

% Initialization

  info.nb_stems            = 0;		% counter of stem vectors
  info.nb_duplicated_stems = 0;		% nb of computed stem vectors that are already in the matrix info.stems
  info.max_stems           = 2*p;	% is increased below when appropriate
  info.stems               = zeros(info.max_stems,p);	% stem vector table in the form of a matrix

  time = tic;

%-------------------------------------------------------------------------------------------------------------------------------
% Compute all the stem vectors of the matrix V without permutating its columns (version with a binary tree)
%-------------------------------------------------------------------------------------------------------------------------------

    bvec      = zeros(p,1);	% initial node of the binary tree
    bvec(1:3) = 1;		% start with bvec = [1 1 1 0 ...] since at least 3 columns of V must be considered (one must have p>2)
    p_bvec    = 3;		% pointer to the current binary tree node
  
    % In the next procedure there is a single QR factorization (inside isf_sv_all_add_from_indices) to compute a null space. The
    % nullity can be deduced from this one

    while p_bvec > 0					% the loop ends when bvec is reset to all 0 by isf_btree at the end of the loop
      sum_bvec = sum(bvec);
      newstem = false;
      if (2 < sum_bvec) && (sum_bvec < r+2)		% at least 3 and at most r+1 columns of V must be considered
        cols = find(bvec~=0);				% select the nonzero indices in bvec, which are the selected columns of V
        [info,newstem] = isf_sv_all_add_from_indices3(V,cols,info);	% add a stem vector to info.stems
        if newstem && (verb2 >= 5)
          % stem in string
          stem = info.stems(info.nb_stems,:);
          stem_str = repmat(' ',p,1);
          stem_str(stem>0) = '1';
          stem_str(stem<0) = '0';
          % bvec in string
          bvec_str = repmat('0',p_bvec,1);
          bvec_str(bvec>0) = '1';
          % print
          fprintf(fout,'\n  | %0s | (%i) from %s',stem_str,info.nb_stems,bvec_str);
        end
      end
      if newstem; bvec(p_bvec) = 0; end			% this prevents exploring the subtree below the current node
      [bvec,p_bvec] = isf_btree(bvec,p_bvec);		% next node of the binary tree
    end

%-------------------------------------------------------------------------------------------------------------------------------
% Purge info.stems of its duplicates
%-------------------------------------------------------------------------------------------------------------------------------

  info.stems               = unique(info.stems(1:info.nb_stems,:),'rows');
  nb_stems                 = size(info.stems,1);
  info.nb_duplicated_stems = info.nb_stems - nb_stems;
  info.nb_stems            = nb_stems;

%-------------------------------------------------------------------------------------------------------------------------------
% Compute the number of signs (nonzero components) for each stem vector
%-------------------------------------------------------------------------------------------------------------------------------

  info.stem_sizes = sum(abs(info.stems),2);

%-------------------------------------------------------------------------------------------------------------------------------
% Printings
%-------------------------------------------------------------------------------------------------------------------------------

  if verb2 >= 2
    fprintf(fout,'\nNumber of stem vectors     = %i',info.nb_stems);
    fprintf(fout,'\nNumber of duplicated stems = %i',info.nb_duplicated_stems);
    fprintf(fout,'\nComputing time             = %g sec',toc(time));
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
      info.stem_zero_indices{i} = I(find(info.stems(I,perm(p+1-i))==0));
    end

  else

    info.stem_zero_indices = [];

  end

  return
