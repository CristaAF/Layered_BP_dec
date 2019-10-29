
N = 128; % VNs; code length 
k = 64; % CNs, i.e., N-k=m (code dimension, length of msg)
m = N-k; %m=rank H numero de filas en H
R=k/N; % code rate R=k/N
clc;
data = load('H.mat')
H = data.H;
[ valid ] = checkValidH2(H,N,m);
BLER = 0;
EbNO=4;
%numCW=1;

[BLER,Itavg] =runScript6( H , R, EbNO );

BLER
Itavg
