function [TRUNCERR,RES,T] = dif1d_orth_it(N,X)
%SOL=@(x) sin(pi*x)+1;
%MSOL=@(x)-pi^2*sin(pi*x);

SOL=@(x) cos(pi*x)-1;
MSOL=@(x)-pi^2*cos(pi*x);

DIF=1.0;
ALPHA=0.9;
XMIN=0.0;
XMAX=1.0;
%N=20;
%KV's

% Immer feiner werdendes Gitter nach Lehrbuch
%X = zeros(1,N+1);
%for I=1:N+1
%  X(I) = XMIN + (ALPHA^(I-1)-1)/(ALPHA^N-1)*(XMAX-XMIN);
%end
%X=linspace(XMIN,XMAX,N+1);

% N halbieren
%for j=1:1 % jeder Durchlauf entfernt Hälfte der Punkte
%  fprintf('abc');
%
%idx =1;
%X2=0;
%for i=1:N+1
%  if mod(i, 2)==1
%    X2(idx)=X(i);
%    idx=idx+1;
%  end
%end
%X=X2;
%N=length(X)-1;
%end

XC = (X(1:N)+X(2:N+1))/2;
XCR = [XMIN, XC, XMAX];

%%% ANALYTISCHE LÖSUNG
TA = zeros(1, N);
for I=1:N
  TA(I)=SOL(XC(I));
end

%figure(1)
%plot(XC, TA, 'x-');
%title('Analytische Loesung')

%%% FVM Lösung

% Randbedingungen
RBE = SOL(XMAX);
RBW = SOL(XMIN);

% Koeffizienten speichern
AP = zeros(1, N);
AE = zeros(1, N);
AW = zeros(1, N);

for I=1:N
  DXE = XCR(I+2)-XCR(I+1);
  DXW = XCR(I+1)-XCR(I);

  AE(I) = DIF/(DXE);
  AW(I) = DIF/(DXW);

  AP(I) = -AE(I)-AW(I);
end

% Gesamtgleichungssystem aufstellen
A = zeros(N);
b = zeros(N,1);

for I=1:N
  b(I) = MSOL(XC(I))*(X(I+1)-X(I));

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

%figure(2)
%plot(XC, T, 'x-');
%title('Numerische Loesung')

%%% Lösungsfehler berechnen
SERR=0.0;
ERR=zeros(1, N);
for I=1:N
    ERR(I)=T(I)-TA(I);
    SERR=SERR+ERR(I)^2;
end
SERR=sqrt(SERR/(N));

%figure(3)
%plot(XC, ERR, 'x-');
%title('Loesungsfehler')

%fprintf('Summierter Fehler %16.10e N=%g\n', SERR, length(X));


%%%% ORDNUNG  BESTIMMEN
% für verdünntes Gitter
% cos, äquidistand
ERR5=3.3333911793e-02;
ERR10=8.2165441917e-03;
ERR20=2.0468563747e-03;
ERR40=5.1125958817e-04;
ERR80=1.2778649795e-04;


%op=log((ERR5)/(ERR10))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );
%op=log((ERR10)/(ERR20))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );
%op=log((ERR20)/(ERR40))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );
%op=log((ERR40)/(ERR80))/log(2);
%fprintf('Ordnung des Verfahrens %16.10e \n',op  );


%%%% RESIDUUM BERECHNEN
t=zeros(N, 1);
for I=1:N
    t(I)=SOL(XC(I));
end

RES=A*t-b;

%figure(4)
%plot(XC,RES,'x-');

%xlabel('XC')
%ylabel('RES')
%title('Residuum')

%%%% Truncation Error berechnen

% b ohne Randwerte
TERR=zeros(N, 1);
for I=1:N
  b(I) = MSOL(XC(I));
end

% Source term
TERRS=zeros(N,1);
for I=1:N
  % benötigte Werte zwischenspeichern
  XP=XCR(I+1);
  Xe=X(I+1);
  Xw=X(I);

  DX = Xe-Xw;
  fP=b(I);

  % west
  if I==1
    fE=b(I+1);
    XE=XCR(I+2);

    TERRS(I) = 1/(6*DX)*((fE-fP)/(XE-XP)-(fP-RBW)/(XP-Xw))*((Xe-XP)^3-(Xw-XP)^3);
  % east
  elseif I==N
    fW=b(I-1);
    XW=XCR(I);

    TERRS(I) = 1/(6*DX)*((RBE-fP)/(Xe-XP)-(fP-fW)/(XP-XW))*((Xe-XP)^3-(Xw-XP)^3);
  % central
  else
    fE=b(I+1);
    XE=XCR(I+2);
    fW=b(I-1);
    XW=XCR(I);

    TERRS(I) = 1/(6*DX)*((fE-fP)/(XE-XP)-(fP-fW)/(XP-XW))*((Xe-XP)^3-(Xw-XP)^3);
  end
end

% Diffusion Term east
TERRDE=zeros(N,1);
TERRDW=zeros(N,1);
for I=1:N
  % benötigte Werte zwischenspeichern
  XP=XCR(I+1);
  TP=T(I);
  fP=b(I);

  Xe=X(I+1);
  Xw=X(I);

  DX = Xe-Xw;

  % west
  if I==1
    XE=XCR(I+2);
    TE=T(I+1);
    fE=b(I+1);
    Xee=X(I+2);

		XEE=XCR(I+3);
    TEE=T(I+2);

    De = (Xe-XP)/(XE-XP);
    Te = De * TE + (1-De) * TP;

    TERRDE(I) = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-RBW)/(XE-Xw))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-RBW)/(XP-Xw)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

    TERRDW(I) = -1/2*((Te-RBW)/(Xe-Xw) - (TP-RBW)/(XP-Xw));

  % west +1
  elseif I==2
    XE=XCR(I+2);
    TE=T(I+1);
    fE=b(I+1);
    Xee=X(I+2);

		XEE=XCR(I+3);
    TEE=T(I+2);

    XW=XCR(I);
    TW=T(I-1);
    fW=b(I-1);
    Xww=X(I-1);

    TERRDE(I) = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

    TERRDW(I) = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-RBW)/(XP-Xww))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-RBW)/(XW-Xww)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));

  % east -1
  elseif I==N-1
    XE=XCR(I+2);
    TE=T(I+1);
    fE=b(I+1);
    Xee=X(I+2);

    XW=XCR(I);
    TW=T(I-1);
    fW=b(I-1);
    Xww=X(I-1);

		XWW=XCR(I-1);
		TWW=T(I-2);

    TERRDE(I) = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((RBE-TE)/(Xee-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

    TERRDW(I) = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-TWW)/(XP-XWW))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));

  % east
  elseif I==N
    XW=XCR(I);
    TW=T(I-1);
    fW=b(I-1);
    Xww=X(I-1);

		XWW=XCR(I-1);
		TWW=T(I-2);

    Dw = (XP-Xw)/(XP-XW);
    Tw = Dw * TW + (1-Dw)*TP;

    TERRDE(I) = 1/2*((RBE-TP)/(Xe-XP) - (RBE-Tw)/(Xe-Xw));

    TERRDW(I) = 1/(2*(XP-XW))*((RBE-TW)/(Xe-XW)-(TP-TWW)/(XP-XWW))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((RBE-TP)/(Xe-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));

  else
    XW=XCR(I);
    TW=T(I-1);
    fW=b(I-1);
    Xww=X(I-1);

    XE=XCR(I+2);
    TE=T(I+1);
    fE=b(I+1);
    Xee=X(I+2);

		XWW=XCR(I-1);
		TWW=T(I-2);

		XEE=XCR(I+3);
    TEE=T(I+2);

    TERRDE(I) = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

    TERRDW(I) = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-TWW)/(XP-XWW))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));
  end
end

for I=1:N
  TERR(I) = TERRS(I)-TERRDE(I)+TERRDW(I);
end


%hold on;
%plot(XC, TERR, 'rx-');

%figure(5)
%plot(XC, RES-TERR, 'x-');
%title('RES - TE')

TRUNCERR=TERR';
end
