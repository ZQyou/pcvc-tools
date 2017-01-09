% Computing Weighted clustering coefficient (edge-weights) for
% vContact weighted network 
%
% INPUT: CSV file path for vContact weighted network
%        CSV file name for the result 
% OUTPUT: Weighted cluster coefficients of genomes
%
function wC=runWCC(csv1,csv2)

  [VC1,W,VC2] = textread(csv1,'%s %f %s','Delimiter',',');
  VC = [VC1; VC2];
  VC = unique(VC);
  N  = size(VC,1);
  NW = size(W,1);
  adj = zeros(N);
  for i = 1:NW
    iVC1 = lookup(VC,VC1{i});
    iVC2 = lookup(VC,VC2{i});
    adj(iVC1,iVC2) = W(i);
  endfor

  adj = adj' + adj;

  wC = weightedClustCoeff(adj);

  fid = fopen(csv2,'w');
  header = {'name','weighted clustering coef'};
  fprintf(fid,'%s,%s\n',header{:});
  formatSpec = '%s,%f\n';
  for row = 1:N
    data = { VC{row}, wC(row) };
    fprintf(fid,formatSpec,data{:});
  end
  fclose(fid);


