close all;
clc;
clear all;


N=20;
ITER=100;

XHIST = zeros(ITER,N+1);
TEHIST = zeros(ITER,N);

% Anfangsgitter
XI = linspace(0,1,N+1);

for A=1:ITER
  XHIST(A,:) = XI;
  % TE ausrechnen
  TERRI = dif1d_orth_it(N,XI);

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
end

% Suche bestes Ergebnis im Vergleich zu 1
figure(1);
hold on;
X1 = XHIST(1,:);
XC1 = (X1(1:N)+X1(2:N+1))/2;
TE1 = TEHIST(1,:);
plot(XC1,TE1,'kx-');
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

figure(2);
hold on;
plot(linspace(1,ITER,ITER), SERR)
plot(linspace(1,ITER,ITER), ERRMIN*ones(1,ITER),'g-')
plot(linspace(1,ITER,ITER), SERR(1)*ones(1,ITER),'r-')

fprintf('Verbesserung relativ: %4.2f%%\n', abs(ERRMIN-SERR(1))/SERR(1)*100);

