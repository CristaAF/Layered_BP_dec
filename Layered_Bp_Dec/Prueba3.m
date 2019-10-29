clear all;
clc;

H=[ 0 0 1 0 1 0 0 1 0 1 0 0;
    1 0 0 0 1 0 0 0 1 0 1 0;
    0 1 0 0 0 1 1 0 0 0 1 0;
    0 0 1 1 0 0 1 0 0 0 0 1;
    1 0 0 0 0 1 0 1 0 0 0 1;
    0 1 0 1 0 0 0 0 1 1 0 0];
   
N = size(H,2);
k=N-size(H,1);
R=k/N;
nbIter=10;
EbNO=4;   % Bit Energy to noise ratio
EsNO  = EbNO+10*log10(R);   %Signal to noise ratio
sigma_ch=sqrt(0.5*((10^(EsNO/10))^-1));   % noise variance, assuming Es=1
LLR_MAX=20;
num_CWs=1; 
CW=0;
Titer=0;
[nbRows,nbCols] = size(H);
x_bits = [1 1 1 1 1 1 0 0 0 0 0 0]; %CW 
x = 2*x_bits -1;   %BPSK modulation
for z=1:num_CWs;
    %CW=CW+1;
y=x+sigma_ch*randn(1,N); 
Lc=2*y/(sigma_ch^2);
Lc=min(Lc,LLR_MAX);
Lc=max(Lc,-LLR_MAX);
counter=0;
Rmj=zeros(1,N);
    SumRmj=zeros(1,N);
%Lc= [1.4683, 4.2594, 1.1861, 2.9966, 0.8355, 0.4020, -2.1164, -0.0846, -2.5058, -2.7922, -0.8071, -1.6935];
LQj=Lc
for nn=1:nbIter
    counter=counter+1;

for n=1:nbRows

Hsub=H(n,:); %get layers of H
[I,J,S] = find(Hsub);%split Hsub row I and column J indices
  J;
  Lqmj=LQj(J) - Rmj(J); %calculate vnd-->cnd message
  Lqmj=min(Lqmj,LLR_MAX);
  Lqmj=max(Lqmj,-LLR_MAX);
  %calculate magnitude of each message(equals beta)
        t=phi(Lqmj);  %  TO Christa: the sum like in 2) t es beta
        t=min(t,LLR_MAX); %clipp again to avoid Inf
        t2=sum(t)-t; %  TO Christa: to exclude i^th element like in 2) %t2 es la suma de betas
        t2=min(t2,LLR_MAX);
        t3=phi(t2);
        t3=min(t3,LLR_MAX);   %  TO Christa: that is now the second part in 2) %te es phi de la suma de betas
        %calculate the sign of each message(equals alpha)
        sgn_s=sign(Lqmj);   %sign of the incoming messages 
        sgn_s(sgn_s==0)=1  %zeros can be set to 1, since magnitude is 0 anyway
        sgn_row=prod(sgn_s); %sign per node % TO Christa:  that is the first part of equation 2)
        
        R=-t3.*sgn_row.*sgn_s; %CND update
        SumRmj(n,J)=real(R);
        Rmj(J)=real(R);
        
        
        
        LQj(J)=Lqmj + real(R);
        LQj=min(LQj,LLR_MAX);
        LQj=max(LQj,-LLR_MAX);
      
end
       LQj=sum(SumRmj);
    %Step5  - calculate output c (hard decision)
        c=sign(LQj+Lc)/2+0.5
        %a=real(c)
        %b=transpose(a)
      %a = real(c)
    %Step 6  - check whether c*H =0, stop if true
        if mod(c*H',2)==0  % %
            break;  %stop decoder
        end

end
it=counter
    if isequal(c,x_bits)
         ErrFrame(z)=0;
    else 
         ErrFrame(z)=1;
    end
    Titer=Titer+counter;
    end
    
BLER = sum(ErrFrame) / num_CWs
AvgIter=Titer/num_CWs

function retVal = phi(x)
% retVal=-log(tanh(abs(x)/2));
retVal=log((exp(x)+1)./(exp(x)-1));
% retVal =  exp(-0.4527*(x.^0.86)+0.0218);
end