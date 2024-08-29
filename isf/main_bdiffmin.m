  clc		% clear command window
  clf
  clear all
  close all	% close all figures
  format long
  format compact

% Add path

  addpath('src','data');

% Set parameters

  verb = 0;	% verbosity level (0: silent, 1: error messages)

% Data

  [A,a,B,b,x] = data_bdiffmin_dgp_2023

% Set options

% options.fout    = fopen('res9','w');
  options.verb    = verb;
  options.bdiffc  = true;
% options.eqtol   = 1.e-8;	% tolerance on the equalities default is 1.e-8

% Run bdiffmin

  init_time = tic;
  [info] = bdiffmin(A,a,B,b,x,options);
  elapsed_time = toc(init_time);
  fprintf('\nElapsed time: %g sec\n',elapsed_time);

  if info.flag
    fprintf('\n\n### BDIFFMIN fails with flag %i\n\n',info.flag);
  end

% Print the B-differential

  fprintf('\nSign vectors in S');
  fprintf('\n^^^^^^^^^^^^^^^^^\n');
  disp(sortrows([info.s;ones(size(info.s))-info.s]));

  fprintf('\nContents of the B-differential');
  fprintf('\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n');
  for i = 1:length(info.bdiff)
    fprintf('Jacobian %i\n',i)
    disp(info.bdiff{i});
  end

% Print the B-differential complement

  if options.bdiffc

    fprintf('\nSign vectors in S complement');
    fprintf('\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n');
    disp(sortrows([info.sc;ones(size(info.sc))-info.sc]));

    fprintf('\nContents of the B-differential complement');
    fprintf('\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n');
    for i = 1:length(info.bdiffc)
      fprintf('Matrix %i\n',i)
      disp(info.bdiffc{i});
    end

  end
