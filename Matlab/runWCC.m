% Computing Weighted clustering coefficient (edge-weights) for
% vContact weighted network 
%
% INPUT: CSV file path for vContact weighted network
%        CSV file name for the result 
% OUTPUT: Weighted cluster coefficients of genomes
%
function [wC,C]=runWCC(csv1,csv2)

  [VC1,W,VC2] = textread(csv1,'%s %f %s','delimiter',',');
% [VC1,W,VC2] = textread(csv1,'%s %f %s','Delimiter',',');
  VC = [VC1; VC2];
  VC = unique(VC);
  N  = size(VC,1);
  NW = size(W,1);
  adj = zeros(N);
  for i = 1:NW
    iVC1 = find(ismember(VC,VC1{i}));
    iVC2 = find(ismember(VC,VC2{i}));
%   iVC1 = lookup(VC,VC1{i});
%   iVC2 = lookup(VC,VC2{i});
    adj(iVC1,iVC2) = W(i);
  end

  adj  = adj' + adj;
  adjW = adj;
  adj(adj~=0) = 1;

% adjW
% adj

  wC = weightedClustCoeff(adjW);
  [avgC,C] = clustCoeff(adj);

  fid = fopen(csv2,'w');
  header = {'name','weighted clustering coef','clustering coef'};
  fprintf(fid,'%s,%s,%s\n',header{:});
  formatSpec = '%s,%f,%f\n';
  for row = 1:N
    data = { VC{row}, wC(row),  C(row) };
    fprintf(fid,formatSpec,data{:});
  end
  fclose(fid);


