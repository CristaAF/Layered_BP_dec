function [ valid ] = checkValidH2(H,N,m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTITUTE OF TELECOMMUNICATIONS  
% University of Stuttgart 
% www.inue-uni-stuttgart.de 
% author: Ahmed Elkelesh
% date: 07.02.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numVNs = isequal(N,size(H,2)); %verifica que N sea igual que el numero de col de H
numCNs = isequal(m,size(H,1)); %verifica que m sea igual que el numero de filas de H
isRateHalf = isequal ( 1 - ( rank(H)/N ) , 0.5 ); %rd=1-(m/n) verify if design rate equal code rate
isBinary = isequal(unique(H),[0;1]); %binario, solo acepta 1 y 0
noRedundVN = all(sum(H,1)~=0);
noRedundCN = all(sum(H,2)~=0);
valid = numVNs & numCNs & isRateHalf & isBinary & noRedundVN & noRedundCN;

end