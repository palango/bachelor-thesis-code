clc
clear all
close all
SOL=@(x,y) sin(pi*x)*sin(pi*y)+1;
MSOL=@(x,y) -2*pi^2*sin(pi*x)*sin(pi*y);
%SOL=@(x,y) sin(pi/2*x)*cos(pi/2*y);
%MSOL=@(x,y) -1*pi^2/2*sin(pi/2*x)*cos(pi/2*y);
DIF=1.0;
XMIN=0.0;
XMAX=1.0;
YMIN=0.0;
YMAX=1.0;
N=5; % KV's in einer Koordinatenrichtung, macht N^2 KV gesamt
NN=N*N;

X = linspace(XMIN, XMAX, N+1);
Y = linspace(YMIN, YMAX, N+1);
XC = (X(1:N)+X(2:N+1))/2;
YC = (Y(1:N)+Y(2:N+1))/2;
XCR = [XMIN, XC, XMAX];
YCR = [YMIN, YC, YMAX];

%%% ANALYTISCHE LÖSUNG
TA = zeros(N, N);
for I=1:N
  for J=1:N
    TA(I, J)=SOL(XC(I), YC(J));
  end
end

figure(1)
surf(XC, YC, TA);
title('Analytische Loesung')
xlabel('X')
ylabel('Y')

%%% FVM Lösung

% Randbedingungen
for I=1:N RBS(I)=SOL(XC(I), 0); end
for I=1:N RBN(I)=SOL(XC(I), 1); end
for I=1:N RBE(I)=SOL(1, YC(I)); end
for I=1:N RBW(I)=SOL(0, YC(I)); end

% Koeffizienten speichern
AP = zeros(N);
AE = zeros(N);
AN = zeros(N);
AW = zeros(N);
AS = zeros(N);

for I=1:N
  for J=1:N
    DX = X(I+1)-X(I);
    DY = Y(J+1)-Y(J);

    DXE = XCR(I+2)-XCR(I+1);
    DXW = XCR(I+1)-XCR(I);

    DYN = YCR(J+2)-YCR(J+1);
    DYS = YCR(J+1)-YCR(J);

    AE(I,J) = DIF*DY/DXE;
    AW(I,J) = DIF*DY/DXW;
    AN(I,J) = DIF*DX/DYN;
    AS(I,J) = DIF*DX/DYS;

    AP(I,J) = -AE(I,J)-AN(I,J)-AW(I,J)-AS(I,J);
  end
end

% Gesamtgleichungssystem aufstellen
A = zeros(NN);
b = zeros(NN, 1);
for J=1:N
  for I=1:N
    IDX = (J-1)*N + I;

    DX = X(I+1)-X(I);
    DY = Y(J+1)-Y(J);
    b(IDX) = MSOL(XC(I),YC(J))*DX*DY;
  end
end

for J=1:N
  for I=1:N
    IDX = (J-1)*N + I;
    % Hauptdiagonale
    A(IDX, IDX) = AP(I,J);

    % Westliche Nebendiagonale
    if mod(IDX,N)==1
      b(IDX) = b(IDX)-AW(I,J)*RBW(J);
    else
      A(IDX, IDX-1) = AW(I,J);
    end

    % Östliche Nebendiagonale
    if mod(IDX,N)==0
      b(IDX) = b(IDX)-AE(I,J)*RBE(J);
    else
      A(IDX, IDX+1) = AE(I,J);
    end

    % Nördliche Nebendiagonale
    if IDX > NN-N
      b(IDX) = b(IDX)-AN(I,J)*RBN(I);
    else
      A(IDX, IDX+N) = AN(I,J);
    end

    % Südliche Nebendiagonale
    if IDX <= N
      b(IDX) = b(IDX)-AS(I,J)*RBS(I);
    else
      A(IDX, IDX-N) = AS(I,J);
    end
  end
end

t=A\b;
T=reshape(t,N,N);

% mit Randwerten
T2 = zeros(N+2);
T2(2:N+1,2:N+1) = T;

figure(2)
surf(XC, YC, T);
title('Numerische Loesung')
xlabel('X')
ylabel('Y')


%%% Lösungsfehler berechnen
SERR=0.0;
ERR=zeros(N);
for I=1:N
  for J=1:N
    ERR(I,J)=T(I, J)-TA(I, J);
    SERR=SERR+ERR(I,J)^2;
  end
end
SERR=sqrt(SERR/NN);

%figure(3)
%surf(XC, YC, ERR);
%title('Loesungsfehler')

fprintf('Summierter Fehler %16.10e NN=%g\n', SERR, NN);

%%%% ORDNUNG  BESTIMMEN
%ERR5=1.6779195551e-02;
%ERR10=4.1327084831e-03;
%ERR20=1.0293533823e-03;
%ERR40=2.5710023908e-04;
%op=log((ERR5)/(ERR10))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );
%op=log((ERR10)/(ERR20))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );
%op=log((ERR20)/(ERR40))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );



%%%% RESIDUUM BERECHNEN
%s=zeros(N);
%for I=1:N
  %for J=1:N
    %s(I, J)=SOL(XC(I), YC(J));
  %end
%end

%s=reshape(s, NN, 1);

%RES=A*s-b;

%RES2 = reshape(RES, N, N);

%figure(4)
%surf(XC, YC, RES2);

%xlabel('XC')
%ylabel('YC')
%zlabel('RES')
%title('Residuum')

%%%% Truncation Error berechnen

%% b ohne Randwerte
%TE=zeros(N);
%b=zeros(N);
%for I=1:N
  %for J=1:N
    %DX = X(I+1)-X(I);
    %DY = Y(J+1)-Y(J);
    %b(IDX) = MSOL(XC(I),YC(J))*DX*DY;
  %end
%end

%for I=3:N-2
  %for J=3:N-2
    %TE(I,J)=1;
    %DX = X(I+1)-X(I);
    %DY = Y(J+1)-Y(J);
    %TE(I,J) = ((DX*DY)/24 * (b(I+1,J)+b(I-1,J)+b(I,J+1)+b(I,J-1)-4*b(I,J)))... % TE_source
            %+ ((T(I+2,J)-3*T(I+1,J)+3*T(I,J)-T(I-1,J))/(24*DX))... % TE_e
            %- ((T(I+1,J)-3*T(I,J)+3*T(I-1,J)-T(I-2,J))/(24*DX))... % TE_w
            %+ ((T(I,J+2)-3*T(I,J+1)+3*T(I,J)-T(I,J-1))/(24*DY))... % TE_n
            %- ((T(I,J+1)-3*T(I,J)+3*T(I,J-1)-T(I,J-2))/(24*DY)); % TE_s
  %end
%end

%% TE Sonderfälle für Randvolumen
%%I=2;
%%DX = X(I+1)-X(I);
%%TE(I) = (DX/24 * (b(I+1) - 2*b(I) + b(I-1)))... % TE_source
      %%+ ((T(I+2)-3*T(I+1)+3*T(I)-T(I-1))/(24*DX))... % TE_e
      %%- ((T(I+1)-3*T(I)+4*T(I-1)-2*RBW)/(24*DX)); % TE_w
%%TE(I) = TE(I)/DX;

%%I = N-1;
%%DX = X(I+1)-X(I);
%%TE(I) = (DX/24 * (b(I+1) - 2*b(I) + b(I-1)))... % TE_source
      %%+ ((2*RBE-4*T(I+1)+3*T(I)-T(I-1))/(24*DX))... % TE_e
      %%- ((T(I+1)-3*T(I)+3*T(I-1)-T(I-2))/(24*DX)); % TE_w
%%TE(I) = TE(I)/DX;

%figure(5)

%surf(XC, YC, TE);
%title('TE');


%RESTE = RES2-TE;
%figure(6)
%surf(XC(3:N-2), YC(3:N-2), RESTE(3:N-2, 3:N-2));
%title('RES-TE');
