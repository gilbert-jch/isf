%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This main program prints on the standard output the tables 5.2 and 5.3
% of the paper 'On the B-differential of the componentwise minimum of
% two affine vector functions' by J.-P. Dussault, J.Ch. Gilbert and B.
% Plaqevent-Jourdain, submitted to 'Mathematical Programming
% Computation'. The run takes much time (several hours or even days).
% For an excerpt of the tables, you can use 'tables_excerpt' instead. To
% run the program, just enter
%
%    tables
%
% in the Matlab window, in the appropriate directory. The output of this
% code obtained on our computer is given in the file 'tables.txt'.
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

  problem_list = ["rand-4-8-2"; "rand-7-8-4"; "rand-7-9-4"; "rand-7-10-5"; "rand-7-11-4"; "rand-7-12-6"; "rand-7-13-5"; "rand-7-14-7"; "rand-8-15-7"; "rand-9-16-8"; "rand-10-17-9";
                  "srand-8-20-2"; "srand-8-20-4"; "srand-8-20-6"; 
                  "rc-2d-20-4"; "rc-2d-20-5"; "rc-2d-20-6"; "rc-2d-20-7"; "rc-2d-20-8";
        	  "rc-perm-5"; "rc-perm-6"; "rc-perm-7"; "rc-perm-8";
        	  "rc-ratio-20-3-7"; "rc-ratio-20-3-9"; "rc-ratio-20-4-7"; "rc-ratio-20-4-9"; "rc-ratio-20-5-7"; "rc-ratio-20-5-9"; "rc-ratio-20-6-7"; "rc-ratio-20-6-9"; "rc-ratio-20-7-7"; "rc-ratio-20-7-9";
                 ];

  solver_list = ["isf_rc"; "isf_a"; "isf_ab"; "isf_abc"; "isf_abcd1"; "isf_abcd2"; "isf_abcd3"; "isf_ad4"];

% Select the output channels and their level of verbosity

  fout  = 1;	% fopen('res','w');
  verb  = 0;	% channel 1 (0: silent, 1: error messages, 2: standard)

  fout2 = 1;	% fopen('res','w');
  verb2 = 0;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)

% Run ISF for each test problem and each solver in the given lists

fout3 = fopen('result','w');
if fout3 < 0
  fprintf('\n### certify_tables: error in opening fout3\n\n');
  return
end
dline = '-------------------------------------------------------------------------------------';
fprintf(fout3,'\n\n\nTable 5.3 in the paper: comparison on the CPU time');
fprintf(fout3,'\n%s',dline);
fprintf(fout3,'\n|                 | Simulated |   ISF (ABCD2)   |   ISF (ABCD3)   |    ISF (AD4)    |');
fprintf(fout3,'\n| Problem         |    RC     |    Time   Ratio |    Time   Ratio |    Time   Ratio |');
fprintf(fout3,'\n%s',dline);

  dline = '------------------------------------------------------------------------------------------------------------------------------------------------------';
  fprintf('\nTable 5.2 in the paper: comparison on the number of LO solved');
  fprintf('\n%s',dline);
  fprintf('\n|                 | Simulated |    ISF (A)     |    ISF (AB)    |   ISF (ABC)    |  ISF (ABCD1)   |  ISF (ABCD2)   |  ISF (ABCD3)   |   ISF (AD4)    |');
  fprintf('\n| Problem         |    RC     |     LOP  Ratio |     LOP  Ratio |     LOP  Ratio |     LOP  Ratio |     LOP  Ratio |     LOP  Ratio |     LOP  Ratio |');
  fprintf('\n%s',dline);

  np = length(problem_list);	% number of problems
  ns = length(solver_list);	% number of solvers

  elapsed_times = zeros(np,ns);
  ratios        = zeros(np,ns);

  for ip = 1:np

    problem = problem_list(ip);

    fprintf('\n| %s%s |    ',problem,repmat(' ',15-strlength(problem),1));

fprintf(fout3,'\n| %s%s |',problem,repmat(' ',15-strlength(problem),1));
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

    % Run all solvers

    for js = 1:ns

      solver = solver_list(js);

      % specify the solver by setting the ISF options

      options = [];
      options.verb  = verb;	% channel 1 (0: silent, 1: error messages, 2: standard)
      options.verb2 = verb;	% channel 2 (0: silent, 1: sign vector, 2: + directions, 3: + intermediate, 4: + verification)
      options.s     = false;	% (list the feasible s's in info.s) true false
      options.sc    = false;	% (list the infeasible s's in info.sc) true false

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
      end

      % Run the ISF solver
      
      init_time = tic;
      [info] = isf(V,options);
      elapsed_times(ip,js) = toc(init_time);

if js == 1
  fprintf(fout3,'  ');
end
if (js == 6) || (js == 7) || (js == 8)
  fprintf(fout3,' %7.2f',elapsed_times(ip,js));
  ratio         = elapsed_times(ip,1)/elapsed_times(ip,js);
  ratios(ip,js) = ratio;
  fprintf(fout3,'  %6.2f |',ratio);
elseif js == 1
  fprintf(fout3,' %7.2f',elapsed_times(ip,js));
  fprintf(fout3,' |');
end

      if js > 1
        if info.flag == 0
          fprintf(' %6i',info.nb_losolve);
        else
          info.nb_losolve = 0;
          fprintf('  -----');
        end
        if info.nb_losolve > 0
          ratio         = nb_losolve_rc/info.nb_losolve;
          ratios(ip,js) = ratio;
          fprintf(' %6.2f | ',ratio);
        else
          ratios(ip,js) = 0;
          fprintf('  ----- | ');
        end
      else
        fprintf('%6i | ',info.nb_losolve);
        nb_losolve_rc = info.nb_losolve;
      end

    end

  end

  fprintf('\n%s',dline);
  ratio_a     = ratios(ratios(:,2)~=0,2);
  ratio_ab    = ratios(ratios(:,3)~=0,3);
  ratio_abc   = ratios(ratios(:,4)~=0,4);
  ratio_abcd1 = ratios(ratios(:,5)~=0,5);
  ratio_abcd2 = ratios(ratios(:,6)~=0,6);
  ratio_abcd3 = ratios(ratios(:,7)~=0,7);
  fprintf('\n| Mean            |           |         %6.2f |         %6.2f |         %6.2f |         %6.2f |         %6.2f |         %6.2f |          ----- |', ...
    mean(ratio_a),mean(ratio_ab),mean(ratio_abc),mean(ratio_abcd1),mean(ratio_abcd2),mean(ratio_abcd3));
  fprintf('\n| Median          |           |         %6.2f |         %6.2f |         %6.2f |         %6.2f |         %6.2f |         %6.2f |          ----- |', ...
    median(ratio_a),median(ratio_ab),median(ratio_abc),median(ratio_abcd1),median(ratio_abcd2),median(ratio_abcd3));
  fprintf('\n%s',dline);

  % Print table 5.3

  dline = '-------------------------------------------------------------------------------------';
  fprintf('\n\n\nTable 5.3 in the paper: comparison on the CPU time');
  fprintf('\n%s',dline);
  fprintf('\n|                 | Simulated |   ISF (ABCD2)   |   ISF (ABCD3)   |    ISF (AD4)    |');
  fprintf('\n| Problem         |    RC     |    Time   Ratio |    Time   Ratio |    Time   Ratio |');
  fprintf('\n%s',dline);

  ratios = zeros(np,ns);

  for ip = 1:length(problem_list)
    problem = problem_list(ip);
    fprintf('\n| %s%s |',problem,repmat(' ',15-strlength(problem),1));
    for js = [1,6,7,8]
      solver = solver_list(js);
      if js == 1
        fprintf('  ');
      end
      fprintf(' %7.2f',elapsed_times(ip,js));
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
  fprintf('\n| Mean            |           |          %6.2f |          %6.2f |          %6.2f |', ...
    mean(ratios(:,6)),mean(ratios(:,7)),mean(ratios(:,8)));
  fprintf('\n| Median          |           |          %6.2f |          %6.2f |          %6.2f |', ...
    median(ratios(:,6)),median(ratios(:,7)),median(ratios(:,8)));
  fprintf('\n%s\n',dline);

  return
