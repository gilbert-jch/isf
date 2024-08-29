%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main program that is used to test the function 'isf'
% (Incremental Sign Feasibility). It is related to the paper 'On the
% B-differential of the componentwise minimum of two affine vectorial
% functions', by Jean-Pierre Duaasault, Jean Charles Gilbert and
% Baptiste Plaquevent-Jourdain. This paper is available on
% http://hal.archives-ouvertes.fr/hal-03872711.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clc		% clear command window
  clf
  clear all
  close all	% close all figures
  format long
  format compact

% Add path

  addpath('src','data');

% Select the output channels and their level of verbosity
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  fout  = 1;	% fopen('res','w');
  verb  = 2;	% channel 1 (0: silent, 1: error messages, 2: standard)

  fout2 = 1;	% fopen('res','w');
  verb2 = 4;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)

% Choose a solver by setting the string 'solver
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Here are the possibilites.

% ''           (empty string, manual tuning below)
% 'isf_rc'     (simulated Rada and Černý's algorithm)
% 'isf_bf'     (brute force)
% 'isf_a'      (combination of the options A, quoted in the paper)
% 'isf_ab'     (combination of the options AB, quoted in the paper)
% 'isf_abc'    (combination of the options ABC, quoted in the paper)
% 'isf_abcd1'  (combination of the options ABCD1, quoted in the paper)
% 'isf_abcd2'  (combination of the options ABCD2, quoted in the paper)
% 'isf_abcd3'  (combination of the options ABCD3, quoted in the paper)
% 'isf_ad4'    (combination of the options AD4, quoted in the paper)

  solver = 'isf_ad4';

% Choose a test-problem among those that have been tested in the paper
% ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% Just uncomment one of the lines below. Each of these files only
% provides an nxp matrix V. Therefore, instead of uncommenting one of
% lines below, providing one's own matrix is perfectly adequate.

% V = data_rand(4,8,2,fout,verb)
% V = data_rand(7,8,4,fout,verb)
% V = data_rand(7,9,4,fout,verb)
% V = data_rand(7,10,5,fout,verb)
% V = data_rand(7,11,4,fout,verb)
% V = data_rand(7,12,6,fout,verb)
% V = data_rand(7,13,5,fout,verb)
% V = data_rand(7,14,7,fout,verb)
% V = data_rand(8,15,7,fout,verb)
% V = data_rand(9,16,8,fout,verb)
% V = data_rand(10,17,9,fout,verb)

% V = data_srand('data_srand_8_20_2.txt',fout,verb)
% V = data_srand('data_srand_8_20_4.txt',fout,verb)
% V = data_srand('data_srand_8_20_6.txt',fout,verb)

% V = data_cerny_rada('data_degen2d_20_4.txt',fout,verb)
% V = data_cerny_rada('data_degen2d_20_5.txt',fout,verb)
% V = data_cerny_rada('data_degen2d_20_6.txt',fout,verb)
% V = data_cerny_rada('data_degen2d_20_7.txt',fout,verb)
% V = data_cerny_rada('data_degen2d_20_8.txt',fout,verb)

% V = data_cerny_rada('data_perm_5.txt',fout,verb)
% V = data_cerny_rada('data_perm_6.txt',fout,verb)
% V = data_cerny_rada('data_perm_7.txt',fout,verb)
% V = data_cerny_rada('data_perm_8.txt',fout,verb)

% V = data_cerny_rada('data_ratio_20_3_7.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_3_9.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_4_7.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_4_9.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_5_7.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_5_9.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_6_7.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_6_9.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_7_7.txt',fout,verb)
% V = data_cerny_rada('data_ratio_20_7_9.txt',fout,verb)

% THAT'S ALL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Other test problems

  V = data_coxeter(4,fout,verb)		% Coxeter arrangement in Rn, one has |S|=n!

% Run algo

  options.fout    = fout;	% fopen('res9','w');
  options.fout2   = fout2;	% fopen('res','w');
  options.verb    = verb;
  options.verb2   = verb2;

  % Run the chosen algorithm

  if strcmp(solver,'isf_bf')		% run the brute force algorithm

    init_time = tic;
    [ns,info] = isf_bf(V,options);
    info.ns = ns;
    elapsed_time = toc(init_time);
    fprintf('\nElapsed time: %g sec (per sign vector: %g sec)\n',elapsed_time,elapsed_time/(2*info.ns));

    return

  elseif strcmp(solver,'isf_rc')	% simulated Rada and Černý's algorithm

    options.rc2018  = true;

  elseif strcmp(solver,'isf_a')		% run ISF (A)

    options.dvnear0 = false;
    options.bestv   = 0;
    options.sv      = 0;

  elseif strcmp(solver,'isf_ab')	% run ISF (AB)

    options.dvnear0 = true;	% Option B
    options.bestv   = 0;
    options.sv      = 0;

  elseif strcmp(solver,'isf_abc')	% run ISF (ABC)

    options.dvnear0 = true;	% Option B
    options.bestv   = 3;	% Option C
    options.sv      = 0;

  elseif strcmp(solver,'isf_abcd1')	% run ISF (ABC3D1)

    options.dvnear0 = true;	% Option B
    options.bestv   = 3;	% Option C
    options.sv      = 1;	% Option D1

  elseif strcmp(solver,'isf_abcd2')	% run ISF (ABC3D2)

    options.dvnear0 = true;	% Option B
    options.bestv   = 3;	% Option C
    options.sv      = 2;	% Option D2

  elseif strcmp(solver,'isf_abcd3')	% run ISF (ABC3D3)

    options.dvnear0 = true;	% Option B
    options.bestv   = 3;	% Option C
    options.sv      = 3;	% Option D3
    options.withd   = true;

  elseif strcmp(solver,'isf_ad4')	% run ISF (AD4); all default options

    options.verb2   = min(verb2,3);

  else					% run ISF (special tuning)

    fprintf('\n\n### main: not recognized solver\n\n');
    return

  end

  options.s  = true;	% (list the feasible s's in info.s) true (default) or false
  options.sc = false;	% info.sc (list of infeasible s's) true false

% Run the ISF solver

  init_time = tic;
  [info] = isf(V,options);
  elapsed_time = toc(init_time);

  if info.flag ~= 0; return; end

  fprintf('\nElapsed time: %g sec (per sign vector: %g sec)\n',elapsed_time,elapsed_time/(2*info.ns));

  if options.s
    fprintf('\nHalf of the %0i feasible sign vectors are\n(the others can be obtained by complementarity):\n\n',info.ns*2)
    disp(info.s)
  end

  if options.sc > 0
    nsc = length(info.sc);
    if options.sc == 1
      fprintf('\nHalf of the %0i infeasible root sign vectors are',nsc*2)
      fprintf('\n(the other infeasible sign vectors can be obtained by completing')
      fprintf('\n'' '' by 0 or 1 [or by setting options.sc = 2] and by complementarity):\n\n')
    else
      fprintf('\nHalf of the %0i infeasible sign vectors are\n(the others can be obtained by complementarity):\n\n',nsc*2)
    end
    for i = 1:nsc
      infosc = info.sc{i};
      for j = 1:length(infosc)
        if isempty(infosc{j})
          fprintf('      ');
        else
          fprintf('     %0i',infosc{j});
        end
      end
      fprintf('\n');
%       disp(info.sc{i});
    end
  end

  return

% Let us do a series of verifications on the output of the function 'isf'

  if info.flag == 0

    fprintf('\nDiagnosis on return from ISF');
    fprintf('\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^');

    [n,m] = size(V);

    % Check that the information in (s,d) solve the problem

    sVd_correct = true;

    for i = 1:info.ns
      s   = 2*blist(i,:)'-1;
      d   = info.d(i,:)';
      Vd  = V'*d;
      sVd = s.*Vd;
      if any(sVd <= 0)
        if sVd_correct; fprintf('\n> s.*(V''*d) is not > 0 for all s and d'); end
        fprintf('\n\n### For s=blist(%i,:) and d=info.d(%i,:), s.*(V''*d)>0 does not hold',i,i);
        fprintf('\n\n     s     V''*d     s.*(V''*d)');
        fprintf('\n    ^^  ^^^^^^^^^^  ^^^^^^^^^');
        for k = 1:m
          fprintf ('\n    %2i  %10.3e  %9.2e',s(k),Vd(k),sVd(k));
        end
        sVd_correct = false;
      end
      if ~sVd_correct; fprintf('\n'); end
    end

    if sVd_correct; fprintf('\n> s.*(V''*d) > 0 for all s and d'); end

    % Check that there is no doubon

    fprintf('\n> There are %i sign vectors',2*info.ns);
    nlist = -ones(2*info.ns,1);
    for i = 1:info.ns
      nlist(i)         = polyval(flip(blist(i,:)),2);
      nlist(info.ns+i) = polyval(flip(1-blist(i,:)),2);
    end
    nlist_purged = unique(nlist);
    if length(nlist_purged) < 2*info.ns
      fprintf('\n\n### Duplicated sign vectors');
    else
      fprintf('\n> No duplicated sign vectors');
    end

    fprintf('\n\n');

    % Put the results in a file

    ilist = zeros(2*info.ns,1);
    for i = 1:info.ns
      ilist(i) = polyval(flip(blist(i,:)),2);
    end
    for i = 1:info.ns
      ilist(info.ns+i) = polyval(flip(1-blist(i,:)),2);
    end
    resfile = fopen('isf_inc_sgn.txt','w');
    fprintf(resfile,'%i\n',sort(ilist));

  end
