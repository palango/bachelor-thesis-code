clc
clear all
close all

SOL=@(x) sin(pi*x)+1;
MSOL=@(x) pi*cos(pi*x)-pi^2*sin(pi*x);

%SOL=@(x) cos(pi*x)-1;
%MSOL=@(x) -pi*sin(pi*x)-pi^2*cos(pi*x);

DIF=1.0;
KONV=1.0;
ALPHA=0.9;
XMIN=0.0;
XMAX=1.0;

N=50; % KV's

X = zeros(1,N+1);
for I=1:N+1
  X(I) = XMIN + (ALPHA^(I-1)-1)/(ALPHA^N-1)*(XMAX-XMIN);
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

  AE(I) = DIF/DXE;
  AW(I) = DIF/DXW;

  AP(I) = -AE(I)-AW(I);

  DCDSE = (X(I+1)-XCR(I+1))/DXE;
  DCDSW = (XCR(I+1)-X(I))/DXW;

  AEK(I) = KONV * DCDSE;
  AWK(I) = -KONV * DCDSW;

  APK(I) = KONV * ((1-DCDSE)-(1-DCDSW));
  AK(I) = KONV;
end

% Gesamtgleichungssystem aufstellen
A = zeros(N);
b = zeros(N,1);

for I=1:N
  b(I) = MSOL(XC(I))*(X(I+1)-X(I));

  % Hauptdiagonale
  A(I,I) = AP(I) + APK(I);

  % Westliche Nebendiagonale
  if mod(I,N)==1
    b(I) = b(I) - AW(I)*RBW + AK(I)*RBW;
  else
    A(I, I-1) = AW(I) + AWK(I);
  end

  % Östliche Nebendiagonale
  if mod(I,N)==0
    b(I) = b(I) - AE(I)*RBE - AK(I)*RBE;
  else
    A(I, I+1) = AE(I) + AEK(I);
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

fprintf('Summierter Fehler %16.10e N=%g\n', SERR, N);


%%% ORDNUNG  BESTIMMEN
ERR5=2.4182329465e-02;
ERR10=5.9638142082e-03;
ERR20=1.4859620816e-03;
ERR40=3.7118037236e-04;
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
TERR=zeros(1, N);
for I=1:N
  b(I) = MSOL(XC(I));
end

for I=3:N-2
  % benötigte Werte zwischenspeichern
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

  % Diffusion Term east
  TERRE = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

  TERRW = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-TWW)/(XP-XWW))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));

  TERRKE = 1/(2*(Xe-Xw))*((TE-TP)/(XE-XP) - (TP-TW)/(XP-XW))*(Xe-XE)*(Xe-XP)...
      + 1/6*(Xe-XP)/DX*((Xe-XP)^2-(XE-XP)^2)...
      * (1/(XE-XP)*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW)) - 1/(XP-XW)*((TE-TW)/(XE-XW) - (TP-TWW)/(XP-XWW)));

  TERRKW = 1/(2*(Xe-Xw))*((TE-TP)/(XE-XP) - (TP-TW)/(XP-XW))*(Xw-XW)*(Xw-XP)...
      + 1/6*(Xw-XP)/DX*((Xw-XP)^2-(XW-XP)^2)...
      * (1/(XE-XP)*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW)) - 1/(XP-XW)*((TE-TW)/(XE-XW) - (TP-TWW)/(XP-XWW)));

  % einzelne Fehler addieren
  TERR(I) = (TERRS-TERRE+TERRW-TERRKE+TERRKW);
end



hold on;
plot(XC, TERR, 'rx-');

RESTE = RES-TERR';

figure(5)
plot(XC(3:N-2), RESTE(3:N-2), 'x-')
xlabel('XC')
ylabel('RES-TE')
