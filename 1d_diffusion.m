% 1D Diffusion Solver
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

function [xc, phi, dx] = diffusion1d(n, f0, fn, f, alpha)
% Anzahl der KV
%n = 25;
x = linspace(0, 1, n + 1);
%x = doublelog(0, 1, n);
xc = zeros(n, 1);
xc(1:n) = 0.5 * (x(1:n) + x(2:n+1));

% delta x
dx = x(2:n+1) - x(1:n);

% Gleichungssystem für äquidistantes Gitter
A = zeros(n);
a = alpha / dx(1)^2;

A = diag(same(-2*a, n)) + diag(same(a, n-1), 1) + diag(same(a, n-1), -1);

% Quellterme in Mittelpunkten berechnen
b = arrayfun(f, xc);

% Randbedingungen in A
A(1, 1) = -3*a;
A(n, n) = -3*a;
% Randbedingungen in b einfügen
b(1) = b(1) - 2*a*f0;
b(n) = b(n) - 2*a*fn;

phi = A\b;

end;%diffusion1d




function [xc, phi, dx] = diffusion1dlog(n, f0, fn, f, alpha)
x = doublelog(0, 1, n);
xc = zeros(n, 1);
xc(1:n) = 0.5 * (x(1:n) + x(2:n+1));

% Größe der KV
dx = diff(x);

% Abstände zwischen Mittelpunkten der KV
rxc = [0, xc', 1];
dp = rxc(2:n+2)-rxc(1:n+1);

% Gleichungssystem für beliebiges Gitter
A = zeros(n);

% Hauptdiagonale
for i=1:n
  A(i,i) = -alpha * (1/(dp(i+1)*dx(i)) + 1/(dp(i)*dx(i)));
endfor;

% Untere Diagonale
for i=2:n
  A(i,i-1) = alpha/(dp(i)*dx(i));
endfor;

% Obere Diagonale
for i=1:n-1
  A(i,i+1) = alpha/(dp(i+1)*dx(i));
endfor;

% Quellterme in Mittelpunkten berechnen
b = arrayfun(f, xc);
b(1) = b(1) - alpha*f0/(dp(1)*dx(1));
b(n) = b(n) - alpha*fn/(dp(n+1)*dx(n));


phi = A\b;

end;%diffusion1dlog











% Konvergenz prüfen
steps = 11;
alpha = 1e-0;
f = @(x)-alpha*pi^2/4*sin(pi/2*x);
i = 1:steps;
% exakte Lösung
solution = @(x)1+sin(pi/2*x);
mean_error = delta_x = deal(zeros(steps, 1));

for j = i
  n = 2^j;
  [xc, phi, dx] = diffusion1dlog(n, 1, 2, f, alpha);
  mean_error(j) = mean(abs(phi- arrayfun(solution, xc)));
  delta_x(j) = dx(1);
endfor;
n=2^3;
% Plotting
[xc1, phi1, dx1] = diffusion1d(n, 1, 2, f, alpha);
[xc2, phi2, dx2] = diffusion1dlog(n, 1, 2, f, alpha);
clf;
hold on;
%plot([0, xc1', 1], [1,hi1', 2], 'bx-')
plot([0, xc2', 1], [1, phi2', 2], 'rx-')
plot([0, xc2', 1], arrayfun(solution, [0, xc2', 1]), 'gx-')
figure;
hold on;
plot([0, xc1', 1], [1, phi1', 2] - arrayfun(solution, [0, xc1', 1]), 'bx-')
plot([0, xc2', 1], [1, phi2', 2] - arrayfun(solution, [0, xc2', 1]), 'rx-')
figure;
semilogy(i, mean_error, '.-.')
xlabel('Anzahl der KV (2^i)')
ylabel('Fehler')

p = log(mean_error(1:steps-1)./mean_error(2:steps))/log(2)
figure;
plot(p);
