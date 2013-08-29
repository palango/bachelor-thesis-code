if 1

% Parameter
alpha = 1e-0;

% Testfunktionen aus Manufactured Solution

%f = @(x)-alpha*pi^2*sin(pi*x);
%solution = @(x)1+sin(pi*x);
%r0=1;
%rn=1;

solution = @(x)-1+cos( pi*x);
f = @(x)-alpha*pi^2*cos(pi*x);
r0 = 0;
rn = -2;

% Anzahl der KV
n=20;
% Plotting
[xc1, phi1, dx1, A1, b1, TE1] = diffusion1d(n, r0, rn, f, alpha);
%[xc2, phi2, dx2, A2, b2] = diffusion1dlog(n, r0, rn, f, alpha);
clf;
figure;
hold on;
title('Loesung');
plot([0, xc1', 1], [r0, phi1', rn], 'bx-')
%plot([0, xc2', 1], [r0, phi2', rn], 'rx-')
%plot([0, xc2', 1], arrayfun(solution, [0, xc2', 1]), 'gx-')

% Fehler
figure;
title('Fehler');
hold on;
error =  [r0, phi1', rn] - arrayfun(solution, [0, xc1', 1]);
plot([0, xc1', 1], error, 'bx-')
%plot([0, xc2', 1], [r0, phi2', rn] - arrayfun(solution, [0, xc2', 1]), 'rx-')

figure;
phi1_exact = arrayfun(solution, xc1);
residuum = (A1*phi1_exact - b1);
title('Residuum');
hold on;
plot(residuum, 'x-')
plot(TE1*n, 'rx')
end;

% Konvergenzprüfung
if 0
% Konvergenz prüfen
%steps = 11;
%i = 1:steps;
%mean_error = delta_x = deal(zeros(steps, 1));

%for j = i
  %n = 2^j;
  %[xc, phi, dx, A, b, TE] = diffusion1dlog(n, 1, 1, f, alpha);
  %mean_error(j) = mean(abs(phi- arrayfun(solution, xc)));
  %delta_x(j) = dx(1);
%endfor;

%figure;
%title('Fehlerkonvergenz');
%semilogy(i, mean_error, '.-.')
%xlabel('Anzahl der KV (2^i)')
%ylabel('Fehler')

%p = log(mean_error(1:steps-1)./mean_error(2:steps))/log(2)

%[haxes,hline1,hline2] = plotyy(1:11, mean_error, 1:10, p, 'semilogy', 'plot');
%axes(haxes(1))
%ylabel('Fehlerkonvergenz')
%axes(haxes(2))
%ylabel('Beobachteter Exponent')
end
