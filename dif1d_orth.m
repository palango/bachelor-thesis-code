clc
clear all
close all

%SOL=@(x) sin(pi*x)+1;
%MSOL=@(x)-pi^2*sin(pi*x);

SOL=@(x) cos(pi*x)-1;
MSOL=@(x)-pi^2*cos(pi*x);

DIF=1.0;
ALPHA=0.9;
XMIN=0.0;
XMAX=1.0;
N=80;% KV's

% Immer feiner werdendes Gitter nach Lehrbuch
X = zeros(1,N+1);
for I=1:N+1
  X(I) = XMIN + (ALPHA^(I-1)-1)/(ALPHA^N-1)*(XMAX-XMIN);
end

% N halbieren
for j=1:1
  fprintf('abc');

idx =1;
X2=0;
for i=1:N+1
  if mod(i, 2)==1
    X2(idx)=X(i);
    idx=idx+1;
  end
end
X=X2;
N=length(X)-1;
end

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

fprintf('Summierter Fehler %16.10e N=%g\n', SERR, length(X));


%%%% ORDNUNG  BESTIMMEN
ERR5=2.4501722712e-02;
ERR10=6.4258655119e-03;
ERR20=2.0705369759e-03;
ERR40=1.0396550077e-03;
ERR80=7.0445840647e-04;

% ALPHA=0.99, cos(x)-2
ERR5=3.3375433977e-02;
ERR10=8.2565958126e-03;
ERR20=2.0862964006e-03;
ERR40=5.4989893414e-04;
ERR80=1.6421419188e-04;

% für verdünntes Gitter
ERR5=1.8313148130e-01;
ERR10=1.1754515384e-01;
ERR20=3.5245466348e-02;
ERR40=8.8115116752e-03;
ERR80=2.1893954570e-03;


op=log((ERR5)/(ERR10))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR10)/(ERR20))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR20)/(ERR40))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR40)/(ERR80))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );


%%%% RESIDUUM BERECHNEN
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

%%%% Truncation Error berechnen

% b ohne Randwerte
TERR=zeros(1, N);
for I=1:N
  b(I) = MSOL(XC(I));
end

for I=3:N-2
  XEE=XCR(I+3);
  XE=XCR(I+2);
  XP=XCR(I+1);
  XW=XCR(I);
  XWW=XCR(I-1);
  Xe=X(I+1);
  Xee=X(I+2);
  Xw=X(I);
  Xww=X(I-1);

  DX = Xe-Xw;

  fE=b(I+1);
  fP=b(I);
  fW=b(I-1);

  TEE=T(I+2);
  TE=T(I+1);
  TP=T(I);
  TW=T(I-1);
  TWW=T(I-2);

  % Source term
  TERRS = (fE-fW)/(2*(XE-XW))*((Xe-XP)^2-(Xw-XP)^2)...
      + 1/(6*DX)*((fE-fP)/(XE-XP)-(fP-fW)/(XP-XW))*((Xe-XP)^3-(Xw-XP)^3);
  %TERRS2 = (DX/24 * (b(I+1) - 2*b(I) + b(I-1)));

  % Diffusion Term east
  TERRE = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));
  %TERRE2 = ((T(I+2)-3*T(I+1)+3*T(I)-T(I-1))/(24*DX));

  TERRW = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-TWW)/(XP-XWW))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));
  %TERRW2 = ((T(I+1)-3*T(I)+3*T(I-1)-T(I-2))/(24*DX));

  TERRC = (TERRS-TERRE+TERRW);
  %TERRC2 = (TERRS2+TERRE2-TERRW2);
  TERR(I) = TERRC/DX; % Division durch Terme bei Quellterm
end


hold on;
plot(XC, TERR, 'rx-');

figure(5)
plot(XC, RES-TERR', 'x-');
title('RES - TE')
