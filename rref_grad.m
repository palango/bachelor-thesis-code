function [X2,W,SERR] = rref_grad(N,X,T)
XC = (X(1:N)+X(2:N+1))/2;
DX = X(2:N+1)-X(1:N);

RB0=0;
RBN=-2;

% gradient berechnen
ERR=zeros(1,N+1);
ERR(1)=(T(1)-RB0)/(DX(1)/2);
for i=2:N
  ERR(i) = (T(i)-T(i-1))/(XC(i)-XC(i-1));
end
ERR(N+1) = (RBN-T(N))/(DX(N)/2);


ERR=abs(ERR);

SUM=0;
for I=1:N+1
  SUM=SUM+ERR(I)^2;
end
SERR=sqrt(SUM/(N+1));

% Verteilung für Ränder der KV berechnen
WX = ERR.^1;
%WX(1)=W(1);
%WX(2:N) = (W(1:N-1)+W(2:N))/2;
%WX(N+1) = W(N);


for C=1:10
  WN=zeros(1,N+1);
  WN(1) = (2*WX(1) + WX(2))/3;
  for I=2:N
    WN(I) = (4*WX(I)+WX(I-1)+WX(I+1))/6;
  end
  WN(N+1)= (2*WX(N+1)+WX(N))/3;
  W=WN;
end

k=(W(1:N)+W(2:N+1))/2;

%figure(1)
%hold on

%plot(linspace(0,1,N+1), ones(1,N+1),'rx')

for I=1:N-1
  XN(1)=0;
  XN(I+1)=(k(I)*X(I) + k(I+1)*X(I+2))/(k(I)+k(I+1));
  XN(N+1)=1;

end
%	plot(XN, (1+C/10)*ones(1,N+1),'bx')

X2=XN;


end %function,1];
