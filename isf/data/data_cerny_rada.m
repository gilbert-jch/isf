function [V] = data_cerny_rada(filename,fout,verb)

%%
% [V] = data_cerny_rada(filename,fout,verb)
%
% Select a problem from Rada and Cerny, "A new algorithm for enumeration
% of cells of hyperplane arrangements and a comparison with Avis and
% Fukuda's reverse search", SIAM Journal on Discrete Mathematics 32:1
% (2018) 455-473 [http://dx.doi.org/10.1137/15M1027930].

  if verb; fprintf(fout,'\nRada and Cerny''s problem ''%s''\n\n',filename); end

  if exist(filename) ~= 2
    if verb; fprintf(fout,'\nUnrecognized file ''%s''\n\n',filename); end
    V = [];
    return
  end

  V = dlmread(filename);	% V = readmatrix(filename) in R2022b

return
