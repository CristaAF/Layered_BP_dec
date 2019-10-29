
function [c,counter] = LDPC_layer_bp_decoder( H, Lc, nbIter)
LLR_MAX=20; 
Lc=min(Lc,LLR_MAX);
Lc=max(Lc,-LLR_MAX);
[nbRows,nbCols] = size(H); 
N=size(H,2);
SumRmj=zeros(1,N);
counter=0;
Rmj=zeros(1,size(H,2));


LQj=Lc;


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
        sgn_s(sgn_s==0)=1;  %zeros can be set to 1, since magnitude is 0 anyway
        sgn_row=prod(sgn_s); %sign per node % TO Christa:  that is the first part of equation 2)
        
        R=-t3.*sgn_row.*sgn_s; %CND update
        SumRmj(n,J)=real(R);%matrix with all previor Rmj
        Rmj(J)=real(R);
        
        LQj(J)=Lqmj + real(R);%VND update
        LQj=min(LQj,LLR_MAX);
        LQj=max(LQj,-LLR_MAX);
      
end
        LQj=sum(SumRmj);%VND values of last subiteration
    %Step5  - calculate output c (hard decision)
        c=sign(LQj+Lc)/2+0.5;
    
    %Step 6  - check whether c*H =0, stop if true
        if mod(c*H',2)==0  % %
            break;  %stop decoder
        end

end
end

function retVal = phi(x)
% retVal=-log(tanh(abs(x)/2));
retVal=log((exp(x)+1)./(exp(x)-1));
% retVal =  exp(-0.4527*(x.^0.86)+0.0218);
end


