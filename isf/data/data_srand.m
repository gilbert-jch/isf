function [V] = data_srand(filename,fout,verb)

%%
% [V] = data_srand(filename,fout,verb)

  if verb; fprintf(fout,'\nStructured random problem ''%s''\n\n',filename); end

  if exist(filename) ~= 2
    if verb; fprintf(fout,'\nUnrecognized file ''%s''\n\n',filename); end
    V = [];
    return
  end

  V = dlmread(filename);	% V = readmatrix(filename) in R2022b

return
