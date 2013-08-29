% FVM für 1D Diffusion mit orthogonalem Gitter
function [xc, phi, dx, A, b] = diffusion1dlog(n, f0, fn, f, alpha)
% Gitter erzeugen
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
