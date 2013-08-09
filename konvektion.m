% 1D Konvektion Solver
function x = same(a, n)
x = a * ones(n, 1);
end;

function x = loginterval(from, to, n)
x = from + (logspace(0, log10(11), n)-1)./(10/(to-from));
end;

function x = doublelog(from, to, n)
half = (to-from)/2;
x = [loginterval(from, half, n/2+1), fliplr(loginterval(to, half, n/2+1))];
x(n/2+1)=[];
end;

function [xc, phi, dx, A, b, TE] = convection1d(n, f0, fn, f, v, rho)

x = linspace(0, 1, n + 1)
xc = zeros(n, 1);
xc(1:n) = 0.5 * (x(1:n) + x(2:n+1));

dx = x(2:n+1) - x(1:n);

% Gleichungssystem für äquidistantes Gitter
A = zeros(n);
a = rho * v / (2 * dx(1));

A = diag(same(0, n)) + diag(same(a, n-1), 1) + diag(same(-a, n-1), -1);

% Quellterme in Mittelpunkten berechnen
b = arrayfun(f, xc)

% Randbedingungen in A
A(1, 1) = a;
A(n, n) = -a;
% Randbedingungen in b einfügen
b(1) = b(1) + 2*a*f0;
b(n) = b(n) - 2*a*fn;

phi = pinv(A)*b;
TE=0;
end;


n=940;
rho = 1;
v = 1;
solution = @(x)-1+cos( pi*x);
f = @(x)-rho*v*pi*sin( pi*x);

[xc1, phi1, dx1, A1, b1, TE1] = convection1d(n, 0,-2, f, rho, v);
plot(xc1, phi1, '-');




%if 1
%% Lösung plotten
%alpha = 1e-0;
%f = @(x)-alpha*pi^2*sin(pi*x);
%% exakte Lösung
%solution = @(x)1+sin(pi*x);

%n=20;
%% Plotting
%[xc1, phi1, dx1, A1, b1, TE1] = diffusion1d(n, 1, 1, f, alpha);
%[xc2, phi2, dx2, A2, b2] = diffusion1dlog(n, 1, 1, f, alpha);
%clf;
%title('Loesung');
%hold on;
%plot([0, xc1', 1], [1, phi1', 1], 'bx-')
%plot([0, xc2', 1], [1, phi2', 1], 'rx-')
%plot([0, xc2', 1], arrayfun(solution, [0, xc2', 1]), 'gx-')

%% Fehler
%figure;
%title('Fehler');
%hold on;
%plot([0, xc1', 1], [1, phi1', 1] - arrayfun(solution, [0, xc1', 1]), 'bx-')
%plot([0, xc2', 1], [1, phi2', 1] - arrayfun(solution, [0, xc2', 1]), 'rx-')

%figure;
%phi1_exact = arrayfun(solution, xc1);
%residuum = abs(A1*phi1_exact - b1);
%title('Residuum');
%semilogy(xc1, residuum, 'x-')
%end;

%if 0
%% Konvergenz prüfen
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
%%title('Fehlerkonvergenz');
%%semilogy(i, mean_error, '.-.')
%%xlabel('Anzahl der KV (2^i)')
%%ylabel('Fehler')

%p = log(mean_error(1:steps-1)./mean_error(2:steps))/log(2)

%[haxes,hline1,hline2] = plotyy(1:11, mean_error, 1:10, p, 'semilogy', 'plot');
%axes(haxes(1))
%ylabel('Fehlerkonvergenz')
%axes(haxes(2))
%ylabel('Beobachteter Exponent')
%end
