function [BLER,Itavg] = runScript6(H,R, EbNO)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INSTITUTE OF TELECOMMUNICATIONS  
% University of Stuttgart 
% www.inue-uni-stuttgart.de 
% author: Ahmed Elkelesh
% date: 07.02.2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = size(H,2);

if nargin<3; EbNO=10; end   % Bit Energy to noise ratio
EbNO
EsNO  = EbNO+10*log10(R);   %Signal to noise ratio
sigma_ch=sqrt(0.5*((10^(EsNO/10))^-1));   % noise variance, assuming Es=1
Iter=0;
num_CWs=1;   %number of codewords, for fast debugging fix to 1 at a high SNR.

maxIter=100;   %Belief propagation maximum number of iterations
% decoder =
% comm.LDPCDecoder('ParityCheckMatrix',sparse(H),'OutputValue','Whole
% codeword','DecisionMethod','Soft
% decision','MaximumIterationCount',maxIter,'NumIterationsOutputPort',1,'IterationTerminationCondition',
% 'Parity check satisfied');  %matlab built-in decoder initialization

ErrFrame = nan(1,num_CWs);
%parfor i=1:num_CWs % parallel cores
for i=1:num_CWs % single core   
    x_bits = zeros(1,N); % we fix our experiment to all-0 codewords.
    x = 2*x_bits -1    %BPSK modulation
    y=x+sigma_ch*randn(1,N);   %AWGN channel 
    Lc=2*y/(sigma_ch^2); %Channel LLR's 
    
%     [LLRxhat, neededIter] = step(decoder, Lc.'); %matlab built-in decoder
%     c=(-sign(LLRxhat))/2+0.5;
    [c,counter] = LDPC_layer_bp_decoder( H, Lc, maxIter)
    if ~isequal(c, x_bits)
        ErrFrame(i)=1;
    else
        ErrFrame(i)=0;
    end
    Iter=Iter+counter;
%     isEqual_all(i)=isequal(c(1:end-1),c2(1:end-1));
end
BLER = sum(ErrFrame) / num_CWs;
Itavg=Iter/num_CWs;
end