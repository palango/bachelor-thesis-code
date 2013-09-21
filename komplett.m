clc
clear all
close all

SOL=@(x)cos(pi*x)-1.0;
BETA=0.5
DIF=1.0;
KONV=1.0;
BV1=SOL(0);
BV2=SOL(1);
XMIN=0.0;
XMAX=1.0;
N=11;

X=zeros(1,N+1);     %Gittereckkoordinaten
XC=zeros(1,N+1);    %Gitterzentren
B=zeros(1,N+1);     %Rechte Seite
b=zeros(1,N-1);     %Rechte Seite
RES=zeros(1,N-1);   %Residuum
A=zeros(N-1,N-1);   %Systemmatrix
T=zeros(1,N+1);     %Lösungsvektor
TA=zeros(1,N+1);    %Analytische Lösung
TE=zeros(1,N+1);    %Lösungsfehler
AE=zeros(1,N+1);
AW=zeros(1,N+1);
AP=zeros(1,N+1);
TE=zeros(1,N-1);

% Gitter
% Gitterabstand
DX=(XMAX-XMIN)/(N-1);

% Gittereckkoordinaten
X(1)=XMIN;
for i=2:N
    X(i)=X(i-1)+DX;
end
X(N+1)=XMAX;

% Gitterzellzentren 
XC(1)=XMIN;
for i=2:N
    XC(i)=X(i-1)+DX/2;
end
XC(N+1)=XMAX;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ANALYTISCHE LÖSUNG 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N+1
    TA(i)=SOL(XC(i));
end
figure(1)
plot(XC,TA,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)

xlabel('XC')
ylabel('TA')
title('Analytische Loesung cos(pi*x)-1.0')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LÖSUNG DIFFUSIONSGLEICHUNG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% DISKRETISIERUNG IM INNEREN DES GEBIETES
for i=2:N-1
    DXPE=XC(i+1)-XC(i);
    %DXWP=XC(i)-XC(i-1);
    DXC=X(i+1)-X(i);
    AE(i)=DIF/(DXPE*DXC);
    AW(i+1)=DIF/(DXPE*DXC);
end

for i=2:N
    AP(i)=(-AE(i)-AW(i));
end

%%% WEST BOUNDARY
D=DIF/((XC(2)-XC(1))*(X(2)-X(1)));
AP(2)=AP(2)-D;
B(2) =B(2) -D*BV1;

%%% EAST BOUNDARY
D=DIF/((XC(N+1)-XC(N))*(X(N)-X(N-1)));
AP(N)=AP(N)-D;
B(N) =B(N) -D*BV2;

%%% QUELLTERM
for i=2:N
    B(i)=B(i)-DIF*pi^2*cos(pi*XC(i));
end

%%% Systemmatrix belegen
% Hauptdiagonale
for i=2:N
    j=i-1;
    A(j,j)=AP(i);
    b(j,1)=B(i);
end
for i=3:N
    j=i-1;
    A(j,j-1)=AW(i);
end
for i=2:N-1
    j=i-1;
    A(j,j+1)=AE(i);
end

%%% LÖSUNG BESTIMMEN
t=A\b;

for i=2:N
   j=i-1;
   T(i)=t(j);
end
T(1)=BV1;
T(N+1)=BV2;

figure(2)
plot(XC,T,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)

xlabel('XC')
ylabel('T')
title('Berechnete Loesung')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% LÖSUNGFEHLER BERECHNEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SERR=0.0;
for i=1:N+1
    ERR(i)=T(i)-TA(i);
    SERR=SERR+ERR(i)^2;
end
SERR=sqrt(SERR/(N-1));

figure(3)
plot(XC,ERR,'--rs','LineWidth',2,...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','g',...
                'MarkerSize',10)

xlabel('XC')
ylabel('ERR')
title('Loesungsfehler')

fprintf('Summierter Fehler %16.10e \n', SERR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ORDNUNG  BESTIMMEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ERR5=5.2635094386e-02;
ERR10=1.0155137137e-02;
ERR20=2.2682751396e-03;
ERR40=5.3782243837e-04; 
op=log((ERR5)/(ERR10))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR10)/(ERR20))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR20)/(ERR40))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% RESIDUUM BERECHNEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=2:N
    j=i-1;
    t(j)=SOL(XC(i));
end

RES=A*t-b;
% RES(1)=0.0;
% RES(N+1)=0.0;

figure(4)
plot(XC(2:N),RES(:,1),'--rs','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10)

xlabel('XC')
ylabel('RES')
title('Residuum')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TRUNCATION ERROR BERECHNEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=3:N-2-1
  DX = X(i+1)-X(i);
  TE_source = DX/24*(SOL(XC(i+1))-2*SOL(XC(i))+SOL(XC(i-1)));
  TE_e = (t(i+2) - 3*t(i+1) + 3*t(i) - t(i-1)) / (24*DX);
  TE_w = (t(i+1) - 3*t(i) + 3*t(i-1) - t(i-2)) / (24*DX);
  TE(i) = TE_source + TE_e - TE_w;
end

% TE Sonderfälle für Randvolumen
%i=2;
%DX = X(i+1)-X(i);
%TE(i) = (DX/24 * (SOL(XC(i+1)) - 2*SOL(XC(i)) + SOL(XC(i-1))))... % TE_source
      %+ ((t(i+2)-3*t(i+1)+3*t(i)-t(i-1))/(24*DX))... % TE_e
      %- ((t(i+1)-3*t(i)+4*t(i-1)-2*SOL(XMIN))/(24*DX)); % TE_w

%i = N-2;
%TE(i) = (DX/24 * (SOL(XC(i+1)) - 2*SOL(XC(i)) + SOL(XC(i-1))))... % TE_source
      %+ ((2*SOL(XMAX)-4*t(i+1)+3*t(i)-t(i-1))/(24*DX))... % TE_e
      %- ((t(i+1)-3*t(i)+3*t(i-1)-t(i-2))/(24*DX)); % TE_w

figure(5)
plot(TE,'--rs','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',10)
title('Truncation Error')
