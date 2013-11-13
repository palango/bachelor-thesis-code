close all;
clc;
clear all;


N=20;
ITER=50;

XHIST = zeros(ITER,N+1);
TEHIST = zeros(ITER,N);
RESHIST = zeros(ITER,N);
THIST = zeros(ITER,N);

% Anfangsgitter
%XI = linspace(0,1,N+1);

%XMIN=0.0;
%XMAX=1.0;
%ALPHA=0.9;
%XI = zeros(1,N+1);
%for I=1:N+1
  %XI(I) = XMIN + (ALPHA^(I-1)-1)/(ALPHA^N-1)*(XMAX-XMIN);
%end
XI=[0,sort(rand(1,N-1)),1];

for A=1:ITER
  XHIST(A,:) = XI;
  % TE ausrechnen
  [TERRI, RESI, TI] = dif1d_orth_it(N,XI);

  ERR=0;
  for I=1:N
    ERR=ERR+TERRI(I)^2;
  end
  SERR(A)=sqrt(ERR/(N));
%  fprintf('Summierter Fehler %16.10e I=%g\n', SERR(A), A);
    

  % Gitter anpassen
  XI = rref2(N,XI,TERRI);

  % Werte speichern
  TEHIST(A,:) = TERRI;
  RESHIST(A,:) =RESI;
THIST(A,:) = TI;
end

% Suche bestes Ergebnis im Vergleich zu 1
figure(1);
hold on;
X1 = XHIST(1,:);
XC1 = (X1(1:N)+X1(2:N+1))/2;
TE1 = TEHIST(1,:);
RES1 = RESHIST(1,:);
plot(XC1,TE1,'kx-');
plot(XC1,RES1,'gx-');
IDXMIN=1;
ERRMIN=SERR(1);
for A=1:ITER
  if SERR(A) <= ERRMIN
    IDXMIN=A;
    ERRMIN=SERR(A);
  end
end

X2 = XHIST(IDXMIN,:);
XC2 = (X2(1:N)+X2(2:N+1))/2;
TE2 = TEHIST(IDXMIN,:);
plot(XC2,TE2,'bx-');
plot(XC2,RESHIST(IDXMIN,:),'rx-');

figure(2);
hold on;
plot(1:ITER, SERR)
plot(1:ITER, ERRMIN*ones(1,ITER),'g-')
plot(1:ITER, SERR(1)*ones(1,ITER),'r-')
fprintf('Verbesserung relativ: %4.2f%%\n', abs(ERRMIN-SERR(1))/SERR(1)*100);

