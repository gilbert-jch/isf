%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This main program prints on the standard output excerpts of the tables
% 5.2 and 5.3 of the paper 'On the B-differential of the componentwise
% minimum of two affine vector functions' by J.-P. Dussault, J.Ch.
% Gilbert and B. Plaqevent-Jourdain, submitted to 'Mathematical
% Programming Computation'. This table excerpt is motivated by the
% desire to make the run fast, while the complete set of problems and
% solvers can take much more time. To run the program, just enter
%
%    tables_excerpt
%
% in the Matlab window, in the appropriate directory. The output of this
% code obtained on our computer is given in the file
% 'tables_excerpt.txt'.
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

  problem_list = ["rand-4-8-2"; "rand-7-8-4"; "rand-7-9-4"; "rand-7-10-5"; "rand-7-11-4";
                  "rc-2d-20-4"; "rc-2d-20-5";
        	  "rc-perm-5";
        	  "rc-ratio-20-3-7"; "rc-ratio-20-3-9"];

  solver_list = ["isf_rc"; "isf_abcd2"; "isf_abcd3"; "isf_ad4"];

% Select the output channels and their level of verbosity

  fout  = 1;	% fopen('res','w');
  verb  = 0;	% channel 1 (0: silent, 1: error messages, 2: standard)

  fout2 = 1;	% fopen('res','w');
  verb2 = 0;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)

% Run ISF for each test problem and each solver in the given lists

  dline = '-------------------------------------------------------------------------';
  fprintf('\nTable 5.2 in the paper (excerpt): comparison on the number of LO solved');
  fprintf('\n%s',dline);
  fprintf('\n|                 | Simulated | ISF (ABCD2) | ISF (ABCD3) |  ISF (AD4)  |');
  fprintf('\n| Problem         |    RC     |  LOP  Ratio |  LOP  Ratio |  LOP  Ratio |');
  fprintf('\n%s',dline);

  np = length(problem_list);	% number of problems
  ns = length(solver_list);	% number of solvers

  elapsed_times = zeros(np,ns);
  ratios        = zeros(np,ns);

  for ip = 1:np

    problem = problem_list(ip);

    fprintf('\n| %s%s |      ',problem,repmat(' ',15-strlength(problem),1));

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
    elseif strcmp(problem,"rc-2d-20-4")
      V = data_cerny_rada('data_degen2d_20_4.txt',fout,verb);
    elseif strcmp(problem,"rc-2d-20-5")
      V = data_cerny_rada('data_degen2d_20_5.txt',fout,verb);
    elseif strcmp(problem,"rc-perm-5")
      V = data_cerny_rada('data_perm_5.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-3-7")
      V = data_cerny_rada('data_ratio_20_3_7.txt',fout,verb);
    elseif strcmp(problem,"rc-ratio-20-3-9")
      V = data_cerny_rada('data_ratio_20_3_9.txt',fout,verb);
    end

    % Run all solvers

    for js = 1:ns

      solver = solver_list(js);

      % specify the solver by setting the ISF options

      options = [];
      options.verb  = verb;	% channel 1 (0: silent, 1: error messages, 2: standard)
      options.verb2 = verb;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)
      options.s     = false;	% (list the feasible s's in info.s) true false (default but used below)

      if strcmp(solver,"isf_rc")
        options.rc2018  = true;
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
      end

      % Run the ISF solver
      
      init_time = tic;
      [info] = isf(V,options);
      elapsed_times(ip,js) = toc(init_time);

      fprintf('%4i',info.nb_losolve);
      if js > 1
        if info.nb_losolve > 0
          ratio         = nb_losolve_rc/info.nb_losolve;
          ratios(ip,js) = ratio;
          fprintf(' %6.2f | ',ratio);
        else
          fprintf('  ----- | ');
        end
      else
        fprintf(' | ');
        nb_losolve_rc = info.nb_losolve;
      end

    end

  end

  fprintf('\n%s',dline);
  ratio_abcd2 = ratios(ratios(:,2)~=0,2);
  ratio_abcd3 = ratios(ratios(:,3)~=0,3);
  fprintf('\n| Mean            |           |      %6.2f |      %6.2f |        ---- |',mean(ratio_abcd2),mean(ratio_abcd3));
  fprintf('\n| Median          |           |      %6.2f |      %6.2f |        ---- |',median(ratio_abcd2),median(ratio_abcd3));
  fprintf('\n%s',dline);

  % Print table 5.3

  dline = '-------------------------------------------------------------------------------';
  fprintf('\n\n\nTable 5.3 in the paper (excerpt): comparison on the CPU time');
  fprintf('\n%s',dline);
  fprintf('\n|                 | Simulated |  ISF (ABCD2)  |  ISF (ABCD3)  |   ISF (AD4)   |');
  fprintf('\n| Problem         |    RC     |  Time   Ratio |  Time   Ratio |  Time   Ratio |');
  fprintf('\n%s',dline);

  ratios = zeros(np,ns);

  for ip = 1:length(problem_list)
    problem = problem_list(ip);
    fprintf('\n| %s%s |',problem,repmat(' ',15-strlength(problem),1));
    for js = 1:length(solver_list)
      solver = solver_list(js);
      if js == 1
        fprintf('    ');
      end
      fprintf(' %5.2f',elapsed_times(ip,js));
      if js > 1
        ratio         = elapsed_times(ip,1)/elapsed_times(ip,js);
        ratios(ip,js) = ratio;
        fprintf('  %6.2f |',ratio);
      else
        fprintf(' |');
      end
    end
  end

  fprintf('\n%s',dline);
  fprintf('\n| Mean            |           |        %6.2f |        %6.2f |        %6.2f |',mean(ratios(:,2)),mean(ratios(:,3)),mean(ratios(:,4)));
  fprintf('\n| Median          |           |        %6.2f |        %6.2f |        %6.2f |',median(ratios(:,2)),median(ratios(:,3)),median(ratios(:,4)));
  fprintf('\n%s\n',dline);

  return
