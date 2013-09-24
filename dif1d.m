clc
clear all
close all

SOL=@(x) sin(pi*x)+1;
MSOL=@(x)-pi^2*sin(pi*x);
DIF=1.0;
XMIN=0.0;
XMAX=1.0;

N=40; % KV's

X = linspace(XMIN, XMAX, N+1);
XC = (X(1:N)+X(2:N+1))/2;
XCR = [XMIN, XC, XMAX];

%%% ANALYTISCHE LÖSUNG
TA = zeros(1, N);
for I=1:N
  TA(I)=SOL(XC(I));
end

figure(1)
plot(XC, TA, 'x-');
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

T=A\b;

figure(2)
plot(XC, T, 'x-');
title('Numerische Loesung')

%%% Lösungsfehler berechnen
SERR=0.0;
ERR=zeros(1, N);
for I=1:N
    ERR(I)=T(I)-TA(I);
    SERR=SERR+ERR(I)^2;
end
SERR=sqrt(SERR/(N));

figure(3)
plot(XC, ERR, 'x-');
title('Loesungsfehler')

fprintf('Summierter Fehler %16.10e \n', SERR);


%%% ORDNUNG  BESTIMMEN
ERR5=2.3729365914e-02;
ERR10=5.8445323862e-03;
ERR20=1.4557255137e-03;
ERR40=3.6359464498e-04;
op=log((ERR5)/(ERR10))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR10)/(ERR20))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR20)/(ERR40))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );


%%% RESIDUUM BERECHNEN
t=zeros(N, 1);
for I=1:N
    t(I)=SOL(XC(I));
end

RES=A*t-b;

figure(4)
plot(XC,RES,'x-');

xlabel('XC')
ylabel('RES')
title('Residuum')

%%% Truncation Error berechnen

% b ohne Randwerte
TE=zeros(1, N);
for I=1:N
  b(I) = MSOL(XC(I));
end

for I=3:N-2
  DX = X(I+1)-X(I);
  TE(I) = (DX/24 * (b(I+1) - 2*b(I) + b(I-1)))... % TE_source
        + ((T(I+2)-3*T(I+1)+3*T(I)-T(I-1))/(24*DX))... % TE_e
        - ((T(I+1)-3*T(I)+3*T(I-1)-T(I-2))/(24*DX)); % TE_w
end;

% TE Sonderfälle für Randvolumen
I=2;
DX = X(I+1)-X(I);
TE(I) = (DX/24 * (b(I+1) - 2*b(I) + b(I-1)))... % TE_source
      + ((T(I+2)-3*T(I+1)+3*T(I)-T(I-1))/(24*DX))... % TE_e
      - ((T(I+1)-3*T(I)+4*T(I-1)-2*RBW)/(24*DX)); % TE_w

I = N-1;
DX = X(I+1)-X(I);
TE(I) = (DX/24 * (b(I+1) - 2*b(I) + b(I-1)))... % TE_source
      + ((2*RBE-4*T(I+1)+3*T(I)-T(I-1))/(24*DX))... % TE_e
      - ((T(I+1)-3*T(I)+3*T(I-1)-T(I-2))/(24*DX)); % TE_w

hold on;
plot(XC, (N-1).*TE, 'rx-');
