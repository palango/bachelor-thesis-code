clc
clear all
close all
SOL=@(x,y) sin(pi*x)*sin(pi*y);
MSOL=@(x,y) -2*pi^2*sin(pi*x)*sin(pi*y);
DIF=1.0;
XMIN=0.0;
XMAX=1.0;
YMIN=0.0;
YMAX=1.0;
N=20; % KV's in einer Koordinatenrichtung, macht N^2 KV gesamt
NN=N*N;

X = linspace(XMIN, XMAX, N+1);
Y = linspace(YMIN, YMAX, N+1);
XC = (X(1:N)+X(2:N+1))/2;
YC = (Y(1:N)+Y(2:N+1))/2;
XCR = [XMIN, XC, XMAX];
YCR = [YMIN, YC, YMAX];

%%% ANALYTISCHE LÖSUNG
Z = zeros(N, N);
for I=1:N
  for J=1:N
    Z(I, J)=SOL(XC(I), YC(J));
  end
end

figure(1)
surf(XC, YC, Z);
title('Analytische Loesung')

%%% FVM Lösung

% Randbedingungen
RBS = zeros(1, N);
for I=1:N RBN(I)=SOL(XC(I), 1); end
for I=1:N RBE(I)=SOL(1, XC(I)); end
RBW = zeros(1, N);

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

    AE(I,J) = DIF/(DXE*DX);
    AW(I,J) = DIF/(DXW*DX);
    AN(I,J) = DIF/(DYN*DY);
    AS(I,J) = DIF/(DYS*DY);

    AP(I,J) = -AE(I,J)-AN(I,J)-AW(I,J)-AS(I,J);
  end
end

% Gesamtgleichungssystem aufstellen
A = zeros(NN);

for J=1:N
  for I=1:N
    IDX = (J-1)*N + I;

    b(IDX) = MSOL(XC(I),YC(J));

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

t=A\b';
T=reshape(t,N,N);

% mit Randwerten
T2 = zeros(N+2);
T2(2:N+1,2:N+1) = T;

figure(2)
surf(XC, YC, T);
title('Numerische Loesung')
