%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This main program can be used to certify the correctedness of the sign
% vectors computed by the solvers 'isf_rc', 'isf_a', 'isf_ab',
% 'isf_abc', 'isf_abcd1', 'isf_abcd2', 'isf_abcd3' and 'isf_ad4'. The
% output of the 'brute force' algorithm serves as a reference in this
% calculation. This code is documented in the paper 'On the
% B-differential of the componentwise minimum of two affine vector
% functions' by J.-P. Dussault, J.Ch. Gilbert and B. Plaqevent-Jourdain,
% submitted to 'Mathematical Programming Computation'. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clc		% clear command window
  clf
  clear all
  close all	% close all figures
  format long
  format compact

% Add path

  addpath('src','data');

% Problem and solver lists

  problem_list = "rand-7-10-5";

  solver_list = ["isf_rc"; "isf_a"; "isf_ab"; "isf_abc"; "isf_abcd1"; "isf_abcd2"; "isf_abcd3"; "isf_ad4"];

% Select the output channels and their level of verbosity

  fout  = 1;	% fopen('res','w');
  verb  = 0;	% channel 1 (0: silent, 1: error messages, 2: standard)

  fout2 = 1;	% fopen('res','w');
  verb2 = 0;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)

% Run ISF for each test problem and each solver in the given lists

  np = length(problem_list);	% number of problems
  ns = length(solver_list);	% number of solvers

  elapsed_times = zeros(np,ns);
  ratios        = zeros(np,ns);

  for ip = 1:np

    problem = problem_list(ip);

    fprintf('Problem ''%s''\n',problem);
    fprintf('^^^^^^^^^%s^',repmat('^',strlength(problem),1));

    % Select one problem

    if strcmp(problem,"rand-4-8-2")
      V = data_rand(4,8,2,fout,verb);
    elseif strcmp(problem,"rand-7-8-4")
      V = data_rand(7,8,4,fout,verb);
    elseif strcmp(problem,"rand-7-9-4")
      V = data_rand(7,9,4,fout,verb);
    elseif strcmp(problem,"rand-7-10-5")
      V = data_rand(7,10,5,fout,verb);
    elseif strcmp(problem,"rand-7-11-4")
      V = data_rand(7,11,4,fout,verb);
    elseif strcmp(problem,"rand-7-12-6")
      V = data_rand(7,12,6,fout,verb);
    elseif strcmp(problem,"rand-7-13-5")
      V = data_rand(7,13,5,fout,verb);
    elseif strcmp(problem,"rand-7-14-7")
      V = data_rand(7,14,7,fout,verb);
    elseif strcmp(problem,"rand-8-15-7")
      V = data_rand(8,15,7,fout,verb);
    elseif strcmp(problem,"rand-9-16-8")
      V = data_rand(9,16,8,fout,verb);
    elseif strcmp(problem,"rand-10-17-9")
      V = data_rand(10,17,9,fout,verb);
    elseif strcmp(problem,"srand-8-20-2")
      V = data_srand('data_srand_8_20_2.txt',fout,verb);
    elseif strcmp(problem,"srand-8-20-4")
      V = data_srand('data_srand_8_20_4.txt',fout,verb);
    elseif strcmp(problem,"srand-8-20-6")
      V = data_srand('data_srand_8_20_6.txt',fout,verb);
    elseif strcmp(problem,"rc-2d-20-4")
      V = data_cerny_rada('data_degen2d_20_4.txt',fout,verb);
    elseif strcmp(problem,"rc-2d-20-5")
      V = data_cerny_rada('data_degen2d_20_5.txt',fout,verb);
    elseif strcmp(problem,"rc-2d-20-6")
      V = data_cerny_rada('data_degen2d_20_6.txt',fout,verb);
    elseif strcmp(problem,"rc-2d-20-7")
      V = data_cerny_rada('data_degen2d_20_7.txt',fout,verb);
    elseif strcmp(problem,"rc-2d-20-8")
      V = data_cerny_rada('data_degen2d_20_8.txt',fout,verb);
    elseif strcmp(problem,"rc-perm-5")
      V = data_cerny_rada('data_perm_5.txt',fout,verb);
    elseif strcmp(problem,"rc-perm-6")
      V = data_cerny_rada('data_perm_6.txt',fout,verb);
    elseif strcmp(problem,"rc-perm-7")
      V = data_cerny_rada('data_perm_7.txt',fout,verb);
    elseif strcmp(problem,"rc-perm-8")
      V = data_cerny_rada('data_perm_8.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-3-7")
      V = data_cerny_rada('data_ratio_20_3_7.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-3-9")
      V = data_cerny_rada('data_ratio_20_3_9.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-4-7")
      V = data_cerny_rada('data_ratio_20_4_7.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-4-9")
      V = data_cerny_rada('data_ratio_20_4_9.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-5-7")
      V = data_cerny_rada('data_ratio_20_5_7.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-5-9")
      V = data_cerny_rada('data_ratio_20_5_9.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-6-7")
      V = data_cerny_rada('data_ratio_20_6_7.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-6-9")
      V = data_cerny_rada('data_ratio_20_6_9.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-7-7")
      V = data_cerny_rada('data_ratio_20_7_7.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-7-9")
      V = data_cerny_rada('data_ratio_20_7_9.txt',fout,verb);
    end

    % Run isf_bf (brute force algorithm) to get the correct list of sign vectors (up to rounding arrors or error made by the
    % optimization solver 'linprog' used by 'bf'). For more information, see section 5.2.1 in the paper 'On the B-differential
    % of the componentwise minimum of two affine vector functions -- The full report'  by J.-P. Dussault, J.Ch. Gilbert and B.
    % Plaqevent-Jourdain [hal-03872711].

    p = size(V,2);

    options = [];
    options.fout  = fout;
    options.verb  = verb;	% channel 1 (0: silent, 1: error messages, 2: standard)
    options.fout2 = fout2;
    options.verb2 = verb2;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)

    fprintf('\nAlgorithm: ''bf'' (provides the reference sign vectors)');

    init_time = tic;
    [info] = bf(V,options);

    elapsed_time = toc(init_time);
    fprintf('\n- elapsed time           = %g sec',elapsed_time);

    nsv = info.ns;
    fprintf('\n- number of sign vectors = %i',nsv*2);
    s_bf = [info.s;ones(nsv,p)-info.s];	% complete with the complementary binary vectors
    s_bf = sortrows(s_bf);				% sort

    % Run all solvers

    for js = 1:ns

      solver = solver_list(js);

      % specify the solver by setting the ISF options

      options = [];
      options.verb  = 0;	% channel 1 (0: silent, 1: error messages, 2: standard)
      options.verb2 = 0;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)
      options.s     = true;	% (list the feasible s's in info.s) true false (default but used below)

      if strcmp(solver,"isf_rc")
        options.rc2018  = true;
      elseif strcmp(solver,"isf_a")
        options.dvnear0 = false;
        options.bestv   = 0;
        options.sv      = 0;
      elseif strcmp(solver,"isf_ab")
        options.dvnear0 = true;	% Option B
        options.bestv   = 0;
        options.sv      = 0;
      elseif strcmp(solver,"isf_abc")
        options.dvnear0 = true;	% Option B
        options.bestv   = 3;	% Option C
        options.sv      = 0;
      elseif strcmp(solver,"isf_abcd1")
        options.dvnear0 = true;	% Option B
        options.bestv   = 3;	% Option C
        options.sv      = 1;	% Option D1
      elseif strcmp(solver,"isf_abcd2")
        options.dvnear0 = true;	% Option B
        options.bestv   = 3;	% Option C
        options.sv      = 2;	% Option D2
      elseif strcmp(solver,"isf_abcd3")
        options.dvnear0 = true;	% Option B
        options.bestv   = 3;	% Option C
        options.sv      = 3;	% Option D3
        options.withd   = true;
      elseif strcmp(solver,"isf_ad4")
        options.sv      = 3;
        options.withd   = false;
      else
        fprintf('\n\n### Unrecognized solver ''%s''\n\n',solver);
        if js < ns
          continue
        else
          break
        end
      end

      % Run the ISF solver
      
      fprintf('\n\nAlgorithm: ''%s''',solver);

      init_time = tic;
      [info] = isf(V,options);
      elapsed_time = toc(init_time);

      if info.flag ~= 0
        fprintf('\n\n### Failure of ISF(%s) on problem ''%s'' (flag = %i)\n\n',solver,problem_list,info.flag);
	return
      end

      fprintf('\n- elapsed time           = %g sec',elapsed_time);

      nsv = info.ns;
      fprintf('\n- number of sign vectors = %i',nsv*2);
      s = [info.s(1:nsv,:);ones(nsv,p)-info.s(1:nsv,:)];	% complete with the complementary binary vectors
      s = sortrows(s);						% sort

      if all(s == s_bf)
        fprintf('\n- the sign vectors are identical to those computed by algorithm ''bf''');
      else
        fprintf('\n- the sign vectors differ from those computed by algorithm ''bf''\n');
        for i = 1:nsv*2
          if any(s(i,:)~=s_bf(i,:))
            fprintf('\n%3i  ',i);
            fprintf('%i',s(i,:));
            fprintf('\n     ');
            fprintf('%i',s_bf(i,:));
          end
        end
      end

    end

    fprintf('\n\n');

  end

  return
