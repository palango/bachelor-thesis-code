clc
clear all
close all

SOL=@(x) sin(pi*x)+1;
MSOL=@(x)-pi^2*sin(pi*x)
DIF=1.0;
XMIN=0.0;
XMAX=1.0;

N=10; % KV's

X = linspace(XMIN, XMAX, N+1);
XC = (X(1:N)+X(2:N+1))/2;
XCR = [XMIN, XC, XMAX];

%%% ANALYTISCHE LÖSUNG
Z = zeros(1, N);
for I=1:N
  Z(I)=SOL(XC(I));
end

figure(1)
plot(XC, Z, 'x-');
title('Analytische Loesung')

%%% FVM Lösung

% Randbedingungen
RBE = SOL(XMAX);
RBW = SOL(XMIN);

% Koeffizienten speichern
AP = zeros(1, N);
AE = zeros(1, N);
AW = zeros(1, N);

for I=1:N
  DX = X(I+1)-X(I);

  DXE = XCR(I+2)-XCR(I+1);
  DXW = XCR(I+1)-XCR(I);

  AE(I) = DIF/(DXE*DX);
  AW(I) = DIF/(DXW*DX);

  AP(I) = -AE(I)-AW(I);
end

% Gesamtgleichungssystem aufstellen
A = zeros(N);
b = zeros(N,1);

for I=1:N
  b(I) = MSOL(XC(I));

  % Hauptdiagonale
  A(I,I) = AP(I);

  % Westliche Nebendiagonale
  if mod(I,N)==1
    b(I) = b(I)-AW(I)*RBW;
  else
    A(I, I-1) = AW(I);
  end

  % Östliche Nebendiagonale
  if mod(I,N)==0
    b(I) = b(I)-AE(I)*RBE;
  else
    A(I, I+1) = AE(I);
  end
end

t=A\b;

figure(2)
plot(XC, t);
title('Numerische Loesung')
