clc
clear all
close all
%SOL=@(x,y) sin(pi*x)*sin(pi*y);
%MSOL=@(x,y) -2*pi^2*sin(pi*x)*sin(pi*y);
SOL=@(x,y) sin(pi/2*x)*sin(pi/2*y);
MSOL=@(x,y) -1*pi^2/2*sin(pi/2*x)*sin(pi/2*y);
DIF=1.0;
XMIN=0.0;
XMAX=1.0;
YMIN=0.0;
YMAX=1.0;
ALPHAX=0.9;
ALPHAY=0.9;
N=40; % KV's in einer Koordinatenrichtung, macht N^2 KV gesamt
NN=N*N;


X = zeros(1,N+1);
for I=1:N+1
  X(I) = XMIN + (ALPHAX^(I-1)-1)/(ALPHAX^N-1)*(XMAX-XMIN);
end

Y = zeros(1,N+1);
for I=1:N+1
  Y(I) = YMIN + (ALPHAY^(I-1)-1)/(ALPHAY^N-1)*(YMAX-YMIN);
end

%NDIVS = 3;
%% N halbieren
%for j=1:NDIVS
%idx =1;
%X2=0;
%for i=1:length(X)+1
  %if mod(i, 2)==1
    %X2(idx)=X(i);
    %idx=idx+1;
  %end
%end
%X=X2;
%N=length(X)-1;
%end

%for j=1:NDIVS
%idx =1;
%Y2=0;
%for i=1:length(Y)+1
  %if mod(i, 2)==1
    %Y2(idx)=Y(i);
    %idx=idx+1;
  %end
%end
%Y=Y2;
%N=length(Y)-1;
%end
%NN=N*N;

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
surf(XC, YC, TA');
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
surf(XC, YC, T');
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

figure(3)
surf(XC, YC, ERR');
title('Loesungsfehler')

fprintf('Summierter Fehler %16.10e NN=%g\n', SERR, NN);

%%% ORDNUNG  BESTIMMEN
ERR5=6.5277704659e-02;
ERR10=1.0773515741e-02;
ERR20=2.2555829695e-03;
ERR40=5.4229194964e-04;
op=log((ERR5)/(ERR10))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR10)/(ERR20))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR20)/(ERR40))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );



%%% RESIDUUM BERECHNEN
s=zeros(N);
for I=1:N
  for J=1:N
    s(I, J)=SOL(XC(I), YC(J));
  end
end

s=reshape(s, NN, 1);

RES=A*s-b;

RES2 = reshape(RES, N, N);

figure(4)
surf(XC, YC, RES2');

xlabel('XC')
ylabel('YC')
zlabel('RES')
title('Residuum')

%%% Truncation Error berechnen

% b ohne Randwerte
TERR=zeros(N);
b=zeros(N);
for I=1:N
  for J=1:N
    DX = X(I+1)-X(I);
    DY = Y(J+1)-Y(J);
    b(I,J) = MSOL(XC(I),YC(J))*DX*DY;
  end
end

for I=3:N-2
  for J=3:N-2
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

    YNN=YCR(J+3);
    YN=YCR(J+2);
    YP=YCR(J+1);
    YS=YCR(J);
    YSS=YCR(J-1);
    Yn=Y(J+1);
    Ynn=Y(J+2);
    Ys=Y(J);
    Yss=Y(J-1);

    DY = Yn-Ys;

    fE=b(I+1,J);
    fP=b(I,J);
    fW=b(I-1,J);
    fN=b(I,J+1);
    fS=b(I,J-1);

    TEE=T(I+2,J);
    TE=T(I+1,J);
    TP=T(I,J);
    TW=T(I-1,J);
    TWW=T(I-2,J);
    TN=T(I,J+1);
    TNN=T(I,J+2);
    TS=T(I,J-1);
    TSS=T(I,J-2);


    DX = X(I+1)-X(I);
    DY = Y(J+1)-Y(J);

    % Source Term
    TERRSO = 1/2*((fE-fP)/((XE-XP)*DX) - (fP-fW)/((XP-XW)*DX))...
           * ((Xe-XP)^3 - (Xw-XP)^3)/3*DY...
          + 1/2*((fN-fP)/((YN-YP)*DY) - (fP-fS)/((YP-YS)*DY))...
           * ((Yn-YP)^3 - (Ys-YP)^3)/3*DX;

    TERRE = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      * (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      + (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      * 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

    TERRW = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-TWW)/(XP-XWW))...
      * (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      + (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      * 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));

    TERRN = 1/(2*(YN-YP))*((TNN-TP)/(YNN-YP)-(TN-TS)/(YN-YS))...
      * (((YP-Yn)^2-(YN-Yn)^2)/(YN-YP))...
      + (1/(Ynn-Yn)*((TNN-TN)/(YNN-YN)-(TN-TP)/(YN-YP)) - 1/(Yn-Ys)*((TN-TP)/(YN-YP)-(TP-TS)/(YP-YS)))...
      * 1/(6*(YN-YP))*(((YP-Yn)^3-(YN-Yn)^3)/(YN-YP));

    TERRS = 1/(2*(YP-YS))*((TN-TS)/(YN-YS)-(TP-TSS)/(YP-YSS))...
      * (((YS-Ys)^2-(YP-Ys)^2)/(YP-YS))...
      + (1/(Yn-Ys)*((TN-TP)/(YN-YP)-(TP-TS)/(YP-YS)) - 1/(Ys-Yss)*((TP-TS)/(YP-YS)-(TS-TSS)/(YS-YSS)))...
      * 1/(6*(YP-YS))*(((YS-Ys)^3-(YP-Ys)^3)/(YP-YS));

    TERR(I,J) = TERRSO - TERRE + TERRW - TERRN + TERRS;
  end
end

figure(5)

surf(XC, YC, TERR');
title('TE');


RESTE = RES2-TERR;
figure(6)
surf(XC(3:N-2), YC(3:N-2), RESTE(3:N-2, 3:N-2)');
title('RES-TE');
