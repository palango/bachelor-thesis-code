close all;
clc;
clear all;


N=20;
ITER=350;

XHIST = zeros(ITER,N+1);
TEHIST = zeros(ITER,N);
RESHIST = zeros(ITER,N);
THIST = zeros(ITER,N);
SOLERRHIST = zeros(ITER,1);
WHIST = zeros(ITER,N+1);
ERRHIST = zeros(ITER,N);


XMIN=0.0;
XMAX=1.0;
ALPHA=0.5;
XI = zeros(1,N+1);
for I=1:N+1
  XI(I) = XMIN + (ALPHA^(I-1)-1)/(ALPHA^N-1)*(XMAX-XMIN);
end
%XI=[0,sort(rand(1,N-1)),1];

for A=1:ITER
  XHIST(A,:) = XI;
  % TE ausrechnen
  [TERRI, RESI, TI,SERRI,ERRI] = dif1d_orth_it(N,XI);


  % Gitter anpassen
  [XI,WI,GRADERRI] = rref_grad(N,XI,TI);
  SERR(A)=GRADERRI;

  % Werte speichern
  TEHIST(A,:) = TERRI;
  RESHIST(A,:) =RESI;
  THIST(A,:) = TI;
  SOLERRHIST(A,:) = SERRI;
  WHIST(A,:) = WI;
  ERRHIST(A,:) = ERRI;
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
IDXMIN=ITER;

X2 = XHIST(IDXMIN,:);
XC2 = (X2(1:N)+X2(2:N+1))/2;
TE2 = TEHIST(IDXMIN,:);
plot(XC2,RESHIST(IDXMIN,:),'rx-');
plot(XC2,TE2,'bx-');
legend('TE1','RES1','TE2','RES2')

figure(2);
hold on;
semilogy(1:ITER, SERR)
semilogy(1:ITER, ERRMIN*ones(1,ITER),'g-')
semilogy(1:ITER, SERR(1)*ones(1,ITER),'r-')
legend('Summierter Gradient','Anfangsfehler','Kleinster Fehler');
fprintf('Verbesserung relativ: %4.2f%%\n', abs(ERRMIN-SERR(1))/SERR(1)*100);

figure(3);
hold on;
semilogy(1:ITER, SOLERRHIST)
semilogy(1:ITER, SOLERRHIST(1)*ones(1,ITER),'g-')
semilogy(1:ITER, SOLERRHIST(IDXMIN)*ones(1,ITER),'r-')
legend('Summierter Fehler','Anfangsfehler','Kleinster Fehler');
fprintf('Verbesserung absoluter Fehler relativ: %4.2f%%\n', abs(SOLERRHIST(1)-SOLERRHIST(IDXMIN))/SOLERRHIST(1)*100);

figure(4);
hold on;
plot(XC1, ERRHIST(1,:),'g-')
plot(XC2, ERRHIST(IDXMIN,:),'r-')
