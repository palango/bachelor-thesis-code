clc
clear all
close all
SOL=@(x,y) sin(pi*x)*sin(pi*y)+1;
MSOL=@(x,y) -2*pi^2*sin(pi*x)*sin(pi*y) ;%+ pi*cos(pi*x)*sin(pi*y) + pi*sin(pi*x)*cos(pi*y);

%SOL=@(x,y) sin(pi*x)*sin(pi*y);
%MSOL=@(x,y) -2*pi^2*sin(pi*x)*sin(pi*y);
DIF=1.0;
XMIN=0.0;
XMAX=1.0;
YMIN=0.0;
YMAX=1.0;
ALPHAX1=0.9;
ALPHAX2=0.8;
ALPHAY1=0.8;
ALPHAY2=0.9;
N=20; % KV's in einer Koordinatenrichtung, macht N^2 KV gesamt
NN=N*N;

% Randwerte
X = zeros(N+1);
for I=1:N+1
  X1 = XMIN + (ALPHAX1^(I-1)-1)/(ALPHAX1^N-1)*(XMAX-XMIN);
  X2 = XMIN + (ALPHAX2^(I-1)-1)/(ALPHAX2^N-1)*(XMAX-XMIN);
  X(I,:) = linspace(X1, X2, N+1);
  %X(:,I)=linspace(XMIN,XMAX,N+1);
end

Y = zeros(N+1);
for I=1:N+1
  Y1 = YMIN + (ALPHAY1^(I-1)-1)/(ALPHAY1^N-1)*(YMAX-YMIN);
  Y2 = YMIN + (ALPHAY2^(I-1)-1)/(ALPHAY2^N-1)*(YMAX-YMIN);
  Y(:,I) = linspace(Y1,Y2,N+1);
  %Y(I,:) = linspace(YMIN,YMAX,N+1);
end

NDIVS = 1;
% N halbieren
for j=1:NDIVS
idx =1;
Xneu=0;
Yneu=0;
for i=1:length(X)
  idy =1;
  inc=false;
  for k=1:length(X)
    if mod(i, 2)==1 & mod(k,2)==1
      Xneu(idx,idy)=X(i,k);
      Yneu(idx,idy)=Y(i,k);
      idy=idy+1;
      inc=true;
    end
  end
  if inc
    idx=idx+1;
  end
end
X=Xneu;
Y=Yneu;
N=length(X)-1;
NN=N*N;
end



% plot mesh
figure(1);
hold on;
grid on;
xlabel('i');
ylabel('j');
for I=1:N+1
  plot(X(I,:), Y(I,:));
  plot(X(:,I), Y(:,I));
end

% cell middle points
XM=zeros(N);
for I=1:N
  for J=1:N
    XM(I,J) = (X(I,J)+X(I+1,J)+X(I,J+1)+X(I+1,J+1))/4;
  end
end
YM=zeros(N);
for I=1:N
  for J=1:N
    YM(I,J) = (Y(I,J)+Y(I+1,J)+Y(I,J+1)+Y(I+1,J+1))/4;
  end
end

%for I=1:N
  %plot(XM(I,:), YM(I,:),'rx');
  %plot(XM(:,I), YM(:,I),'rx');
%end

XMR=zeros(N+2);
XMR(2:N+1,2:N+1) = XM;
XMR(1,:) = XMIN*ones(1,N+2);
XMR(N+2,:) = XMAX*ones(1,N+2);
for I=1:N
    XMR(I+1,1) = (X(I,1)+X(I+1,1))/2;
end
for I=1:N
    XMR(I+1,N+2) = (X(I,N+1)+X(I+1,N+1))/2;
end

YMR=zeros(N+2);
YMR(2:N+1,2:N+1) = YM;
YMR(:,1) = YMIN*ones(N+2,1);
YMR(:,N+2) = YMAX*ones(N+2,1);
for J=1:N
    YMR(1,J+1) = (Y(1,J)+Y(1,J+1))/2;
end
for J=1:N
    YMR(N+2,J+1) = (Y(N+1,J)+Y(N+1,J+1))/2;
end

for I=1:N+2
  plot(XMR(I,:), YMR(I,:),'rx');
  plot(XMR(:,I), YMR(:,I),'rx');
end

% metrics
MXXIH =zeros(N+1,N);
MXETAH=zeros(N+1,N);
MYXIH =zeros(N+1,N);
MYETAH=zeros(N+1,N);

MXXIV =zeros(N,N+1);
MXETAV=zeros(N,N+1);
MYXIV =zeros(N,N+1);
MYETAV=zeros(N,N+1);

for I=1:N+1
  for J=1:N
    X1=XMR(I,J+1);
    X2=XMR(I+1,J+1);
    Y1=YMR(I,J+1);
    Y2=YMR(I+1,J+1);

    XN=X(I,J+1);
    YN=Y(I,J+1);
    XS=X(I,J);
    YS=Y(I,J);

    NORM=sqrt((X2-X1)^2 + (Y2-Y1)^2);
    MXXIH(I,J)=(X2-X1)/NORM;
    MYXIH(I,J)=(Y2-Y1)/NORM;

    DS=sqrt((XN-XS)^2+(YN-YS)^2);
    MXETAH(I,J)=(XN-XS)/DS;
    MYETAH(I,J)=(YN-YS)/DS;
  end
end
for I=1:N
  for J=1:N+1
    X1=XMR(I+1,J);
    X2=XMR(I+1,J+1);
    Y1=YMR(I+1,J);
    Y2=YMR(I+1,J+1);

    XE=X(I+1,J);
    YE=Y(I+1,J);
    XW=X(I,J);
    YW=Y(I,J);

    NORM=sqrt((X2-X1)^2 + (Y2-Y1)^2);
    MXETAV(I,J)=(X2-X1)/NORM;
    MYETAV(I,J)=(Y2-Y1)/NORM;

    DS=sqrt((XE-XW)^2+(YE-YW)^2);
    MXXIV(I,J)=(XE-XW)/DS;
    MYXIV(I,J)=(YE-YW)/DS;
  end
end

% cell face centers
% north-south sides
XCV=zeros(N,N+1);
YCV=zeros(N,N+1);
YCV(:,N+1)=Y(1:N,N+1);
for I=1:N
  for J=1:N+1
    XCV(I,J) = (X(I,J)+X(I+1,J))/2;
    YCV(I,J) = (Y(I,J)+Y(I+1,J))/2;
  end
end
for I=1:N+1
  plot(XCV(:,I), YCV(:,I),'kx');
end

% east-west sides
XCH=zeros(N+1,N);
XCH(N+1,:)=X(N+1,1:N);
YCH=zeros(N+1,N);
for I=1:N+1
  for J=1:N
    XCH(I,J) = (X(I,J)+X(I,J+1))/2;
    YCH(I,J) = (Y(I,J)+Y(I,J+1))/2;
  end
end
for I=1:N+1
  plot(XCH(I,:), YCH(I,:),'gx');
end

%% normal vectors
DXH = zeros(N+1,N);
DYH = zeros(N+1, N);
LENGTHH = zeros(N+1, N);
NVHX = zeros(N+1,N);
NVHY = zeros(N+1,N);
for I=1:N+1
  for J=1:N
    DX = X(I,J+1)-X(I,J);
    DY = Y(I,J+1)-Y(I,J);

    NORM = sqrt(DX^2 + DY^2);
    NVHX(I,J) = DY/NORM;
    NVHY(I,J) = -DX/NORM;

    DXH(I,J)=DX;
    DYH(I,J)=DY;
    LENGTHH(I,J)=NORM;
  end
end
quiver(XCH, YCH, NVHX, NVHY, 0.2);

DXV = zeros(N,N+1);
DYV = zeros(N,N+1);
LENGTHV = zeros(N,N+1);
NVVX = zeros(N,N+1);
NVVY = zeros(N,N+1);
for I=1:N
  for J=1:N+1
    DX = X(I+1,J)-X(I,J);
    DY = Y(I+1,J)-Y(I,J);

    NORM = sqrt(DX^2 + DY^2);
    NVVX(I,J) = -DY/NORM;
    NVVY(I,J) = DX/NORM;

    DXV(I,J)=DX;
    DYV(I,J)=DY;
    LENGTHV(I,J)=NORM;
  end
end
quiver(XCV, YCV, NVVX, NVVY, 0.2);

DMH=zeros(N+1,J);
MJH=zeros(N+1,N);
for I=1:N+1
  for J=1:N
    X1=XMR(I,J+1);
    X2=XMR(I+1,J+1);
    Y1=YMR(I,J+1);
    Y2=YMR(I+1,J+1);

    NORM=sqrt((X2-X1)^2 + (Y2-Y1)^2);
    DMH(I,J)=NORM;

    XN=X(I,J+1);
    YN=Y(I,J+1);
    XS=X(I,J);
    YS=Y(I,J);

    TOP = (X2-X1)*(YN-YS) - (Y2-Y1)*(XN-XS);

    MJH(I,J) = TOP / (NORM * LENGTHH(I,J));
  end
end
DMV=zeros(N,N+1);
MJV=zeros(N,N+1);
for I=1:N
  for J=1:N+1
    X1=XMR(I+1,J);
    X2=XMR(I+1,J+1);
    Y1=YMR(I+1,J);
    Y2=YMR(I+1,J+1);

    NORM=sqrt((X2-X1)^2 + (Y2-Y1)^2);
    DMV(I,J)=NORM;

    XE=X(I+1,J);
    YE=Y(I+1,J);
    XW=X(I,J);
    YW=Y(I,J);

    TOP=(XE-XW)*(Y2-Y1) - (YE-YW)*(X2-X1);

    MJV(I,J) = TOP / (NORM * LENGTHV(I,J));
  end
end

% cell volume
V=zeros(N);
for I=1:N
  for J=1:N
    Xne = X(I+1,J+1);
    Xnw = X(I,J+1);
    Xsw = X(I,J);
    Xse = X(I+1,J);

    Yne = Y(I+1,J+1);
    Ynw = Y(I,J+1);
    Ysw = Y(I,J);
    Yse = Y(I+1,J);

    V(I,J) = 0.5*abs((Xse-Xnw)*(Yne-Ysw) - (Xne-Xsw)*(Yse-Ynw));
    text(X(I,J)+0.01,Ysw+0.01,num2str(V(I,J)));
  end
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

%XC = (X(1:N)+X(2:N+1))/2;
%YC = (Y(1:N)+Y(2:N+1))/2;
%XCR = [XMIN, XC, XMAX];
%YCR = [YMIN, YC, YMAX];

%%%% ANALYTISCHE LÖSUNG
TA = zeros(N, N);
for I=1:N
  for J=1:N
    TA(I, J)=SOL(XM(I,J), YM(I,J));
  end
end

hold off;
figure(9)
surf(XM, YM, TA');
title('Analytische Loesung')
xlabel('X')
ylabel('Y')

%%%% FVM Lösung

%% Randbedingungen
RB=zeros(N+2);
for I=1:N+2
  for J=1:N+2
    RB(I,J) = SOL(XMR(I,J),YMR(I,J));
  end
end
RB(2:N+1,2:N+1) = zeros(N);

%%% Koeffizienten speichern
AP = zeros(N);
AE = zeros(N);
AN = zeros(N);
AW = zeros(N);
AS = zeros(N);
ANE= zeros(N);
ANW= zeros(N);
ASE= zeros(N);
ASW= zeros(N);
for I=1:N
  for J=1:N
    % normal terms, orth grid
    PDE = MYETAH(I+1,J)*NVHX(I+1,J) - MXETAH(I+1,J)*NVHY(I+1,J);
    PDW = MYETAH(I,J)*NVHX(I,J) - MXETAH(I,J)*NVHY(I,J);

    PDN = MXXIV(I,J+1)*NVVY(I,J+1) - MYXIV(I,J+1)*NVVX(I,J+1);
    PDS = MXXIV(I,J)*NVVY(I,J) - MYXIV(I,J)*NVVX(I,J);

    AE(I,J) = DIF*LENGTHH(I+1,J) * PDE/MJH(I+1,J)/DMH(I+1,J);
    AW(I,J) = DIF*LENGTHH(I,J) * PDW/MJH(I,J)/DMH(I,J);
    AN(I,J) = DIF*LENGTHV(I,J+1) * PDN/MJV(I,J+1)/DMV(I,J+1);
    AS(I,J) = DIF*LENGTHV(I,J) * PDS/MJV(I,J)/DMV(I,J);

    AP(I,J) = -AE(I,J)-AN(I,J)-AW(I,J)-AS(I,J);

    % terms from unorthogonal grid
    PDE2 = MXXIH(I+1,J)*NVHY(I+1,J) - MYXIH(I+1,J)*NVHX(I+1,J);
    PDW2 = MXXIH(I,J)*NVHY(I,J) - MYXIH(I,J)*NVHX(I,J);

    PDN2 = MYETAV(I,J+1)*NVVX(I,J+1) - MXETAV(I,J+1)*NVVY(I,J+1);
    PDS2 = MYETAV(I,J)*NVVX(I,J) - MXETAV(I,J)*NVVY(I,J);

    % east
    PDE3 = (DIF*LENGTHH(I+1,J)*PDE2/MJH(I+1,J)/LENGTHH(I+1,J))/4;
    AN(I,J) = AN(I,J) + PDE3;
    ANE(I,J) = PDE3;
    AS(I,J) = AS(I,J) - PDE3;
    ASE(I,J) = -PDE3;

    % west
    PDW3 = (DIF*LENGTHH(I,J)*PDW2/MJH(I,J)/LENGTHH(I,J))/4;
    AN(I,J) = AN(I,J) + PDW3;
    ANW(I,J) = PDW3;
    AS(I,J) = AS(I,J) - PDW3;
    ASW(I,J) = -PDW3;

    % north
    PDN3 = (DIF*LENGTHV(I,J+1)*PDN2/MJV(I,J+1)/LENGTHV(I,J+1))/4;
    AE(I,J)  = AE(I,J) + PDN3;
    ANE(I,J) = ANE(I,J)+ PDN3;
    AW(I,J)  = AW(I,J) - PDN3;
    ANW(I,J) = ANW(I,J)- PDN3;

    % south
    PDS3 = (DIF*LENGTHV(I,J)*PDS2/MJV(I,J)/LENGTHV(I,J))/4;
    AE(I,J)  = AE(I,J) + PDS3;
    ASE(I,J) = ASE(I,J)+ PDS3;
    AW(I,J)  = AW(I,J) - PDS3;
    ASW(I,J) = ASW(I,J)- PDS3;
  end
end

%% Gesamtgleichungssystem aufstellen
A = zeros(NN);
b = zeros(NN, 1);
for J=1:N
  for I=1:N
    IDX = (J-1)*N + I;
    b(IDX) = MSOL(XM(I,J),YM(I,J))*V(I,J);
  end
end

for J=1:N
  for I=1:N
    IDX = (J-1)*N + I;

    % Hauptdiagonale
    A(IDX, IDX) = AP(I,J);

    % Westliche Nebendiagonale
    if mod(IDX,N)==1
      b(IDX) = b(IDX)-AW(I,J)*RB(1,J+1);
    else
      A(IDX, IDX-1) = AW(I,J);
    end

    % Östliche Nebendiagonale
    if mod(IDX,N)==0
      b(IDX) = b(IDX)-AE(I,J)*RB(N+2,J+1);
    else
      A(IDX, IDX+1) = AE(I,J);
    end

    % Nördliche Nebendiagonale
    if IDX > NN-N
      b(IDX) = b(IDX)-AN(I,J)*RB(I+1,N+2);
    else
      A(IDX, IDX+N) = AN(I,J);
    end

    % Südliche Nebendiagonale
    if IDX <= N
      b(IDX) = b(IDX)-AS(I,J)*RB(I+1,1);
    else
      A(IDX, IDX-N) = AS(I,J);
    end

    % non-orthogonal factors
    % south-west
    %if mod(IDX,N)==1 | IDX <= N
      %b(IDX) = b(IDX)-ASW(I,J)*RB(I,J);
    %else
      %A(IDX,IDX-N-1) = ASW(I,J);
    %end

    %% north-west
    %if mod(IDX,N)==1 | IDX > NN-N
      %b(IDX) = b(IDX)-ANW(I,J)*RB(I,J+2);
    %else
      %A(IDX,IDX+N-1) = ANW(I,J);
    %end

    %% south-east
    %if mod(IDX,N)==0 | IDX <= N
      %b(IDX) = b(IDX)-ASE(I,J)*RB(I+2,J);
    %else
      %A(IDX,IDX-N+1) = ASE(I,J);
    %end

    %% north-east
    %if mod(IDX,N)==0 | IDX > NN-N
      %b(IDX) = b(IDX)-ANE(I,J)*RB(I+2,J+2);
    %else
      %A(IDX,IDX+N+1) = ANE(I,J);
    %end
  end
end

t=A\b;
T=reshape(t,N,N);

% mit Randwerten
T2 = RB;
T2(2:N+1,2:N+1) = T;

figure(2)
surf(XM, YM, T');
title('Numerische Loesung1');
xlabel('X')
ylabel('Y')


%%%% Lösungsfehler berechnen
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
surf(XM, YM, ERR');
title('Loesungsfehler')

fprintf('Summierter Fehler %16.10e NN=%g\n', SERR, NN);

%%%% ORDNUNG  BESTIMMEN
ERR5=5.9670806291e-02;
ERR10=6.1601523373e-02;
ERR20=9.8536823766e-02;
%ERR40=5.4229194964e-04;
op=log((ERR5)/(ERR10))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
op=log((ERR10)/(ERR20))/log(2);
fprintf('Ordnung des Verfahrens %16.10e \n',op  );
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
%surf(XC, YC, RES2');

%xlabel('XC')
%ylabel('YC')
%zlabel('RES')
%title('Residuum')

%%%% Truncation Error berechnen

%% b ohne Randwerte
%TERR=zeros(N);
%b=zeros(N);
%for I=1:N
  %for J=1:N
    %DX = X(I+1)-X(I);
    %DY = Y(J+1)-Y(J);
    %b(I,J) = MSOL(XC(I),YC(J))*DX*DY;
  %end
%end

%for I=3:N-2
  %for J=3:N-2
    %% benötigte Werte zwischenspeichern
    %XEE=XCR(I+3);
    %XE=XCR(I+2);
    %XP=XCR(I+1);
    %XW=XCR(I);
    %XWW=XCR(I-1);
    %Xe=X(I+1);
    %Xee=X(I+2);
    %Xw=X(I);
    %Xww=X(I-1);

    %DX = Xe-Xw;

    %YNN=YCR(J+3);
    %YN=YCR(J+2);
    %YP=YCR(J+1);
    %YS=YCR(J);
    %YSS=YCR(J-1);
    %Yn=Y(J+1);
    %Ynn=Y(J+2);
    %Ys=Y(J);
    %Yss=Y(J-1);

    %DY = Yn-Ys;

    %fE=b(I+1,J);
    %fP=b(I,J);
    %fW=b(I-1,J);
    %fN=b(I,J+1);
    %fS=b(I,J-1);

    %TEE=T(I+2,J);
    %TE=T(I+1,J);
    %TP=T(I,J);
    %TW=T(I-1,J);
    %TWW=T(I-2,J);
    %TN=T(I,J+1);
    %TNN=T(I,J+2);
    %TS=T(I,J-1);
    %TSS=T(I,J-2);


    %DX = X(I+1)-X(I);
    %DY = Y(J+1)-Y(J);

    %% Source Term
    %TERRSO = 1/2*((fE-fP)/((XE-XP)*DX) - (fP-fW)/((XP-XW)*DX))...
           %* ((Xe-XP)^3 - (Xw-XP)^3)/3*DY...
          %+ 1/2*((fN-fP)/((YN-YP)*DY) - (fP-fS)/((YP-YS)*DY))...
           %* ((Yn-YP)^3 - (Ys-YP)^3)/3*DX;

    %TERRE = 1/(2*(XE-XP))*((TEE-TP)/(XEE-XP)-(TE-TW)/(XE-XW))...
      %* (((XP-Xe)^2-(XE-Xe)^2)/(XE-XP))...
      %+ (1/(Xee-Xe)*((TEE-TE)/(XEE-XE)-(TE-TP)/(XE-XP)) - 1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)))...
      %* 1/(6*(XE-XP))*(((XP-Xe)^3-(XE-Xe)^3)/(XE-XP));

    %TERRW = 1/(2*(XP-XW))*((TE-TW)/(XE-XW)-(TP-TWW)/(XP-XWW))...
      %* (((XW-Xw)^2-(XP-Xw)^2)/(XP-XW))...
      %+ (1/(Xe-Xw)*((TE-TP)/(XE-XP)-(TP-TW)/(XP-XW)) - 1/(Xw-Xww)*((TP-TW)/(XP-XW)-(TW-TWW)/(XW-XWW)))...
      %* 1/(6*(XP-XW))*(((XW-Xw)^3-(XP-Xw)^3)/(XP-XW));

    %TERRN = 1/(2*(YN-YP))*((TNN-TP)/(YNN-YP)-(TN-TS)/(YN-YS))...
      %* (((YP-Yn)^2-(YN-Yn)^2)/(YN-YP))...
      %+ (1/(Ynn-Yn)*((TNN-TN)/(YNN-YN)-(TN-TP)/(YN-YP)) - 1/(Yn-Ys)*((TN-TP)/(YN-YP)-(TP-TS)/(YP-YS)))...
      %* 1/(6*(YN-YP))*(((YP-Yn)^3-(YN-Yn)^3)/(YN-YP));

    %TERRS = 1/(2*(YP-YS))*((TN-TS)/(YN-YS)-(TP-TSS)/(YP-YSS))...
      %* (((YS-Ys)^2-(YP-Ys)^2)/(YP-YS))...
      %+ (1/(Yn-Ys)*((TN-TP)/(YN-YP)-(TP-TS)/(YP-YS)) - 1/(Ys-Yss)*((TP-TS)/(YP-YS)-(TS-TSS)/(YS-YSS)))...
      %* 1/(6*(YP-YS))*(((YS-Ys)^3-(YP-Ys)^3)/(YP-YS));

    %TERR(I,J) = TERRSO - TERRE + TERRW - TERRN + TERRS;
  %end
%end

%figure(5)

%surf(XC, YC, TERR');
%title('TE');


%RESTE = RES2-TERR;
%figure(6)
%surf(XC(3:N-2), YC(3:N-2), RESTE(3:N-2, 3:N-2)');
%title('RES-TE');
