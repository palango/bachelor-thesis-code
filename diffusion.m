% 1D Diffusion Solver
1; 

% Vektor von n gleichen Elementen erzeugen
function x = same(a, n)
x = a * ones(n, 1);
end;

% Vector von logarithmisch verteilen Elementen in [from, to] erzeugen
function x = loginterval(from, to, n)
x = from + (logspace(0, log10(11), n)-1)./(10/(to-from));
end;

% Vektor von doppelt log. verteilten Elementen mit höherer Dichte am Rand
function x = doublelog(from, to, n)
half = (to-from)/2;
x = [loginterval(from, half, n/2+1), fliplr(loginterval(to, half, n/2+1))];
x(n/2+1)=[];
end;

% FVM für 1D Diffusion mit äquidistantem Gitter
function [xc, phi, dx, A, b, TE] = diffusion1d(n, f0, fn, f, alpha)

% Gitter erzeugen
x = linspace(0, 1, n + 1);
% Mittelpunkte der KV
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
% Randbedingungen in b einfügen,
% b ohne RB wirde für TE benötigt
b_rb = b;
b_rb(1) = b_rb(1) - 2*a*f0;
b_rb(n) = b_rb(n) - 2*a*fn;

phi = A\b_rb;


% Truncation Error berechnen
for i=3:n-2
  TE(i) = (dx(1)/24 * (b(i+1) - 2*b(i) + b(i-1)))... % TE_source
        + ((phi(i+2)-3*phi(i+1)+3*phi(i)-phi(i-1))/(24*dx(1)))... % TE_e
        - ((phi(i+1)-3*phi(i)+3*phi(i-1)-phi(i-2))/(24*dx(1))); % TE_w
endfor;

% TE Sonderfälle für Randvolumen
i=2;
TE(i) = (dx(1)/24 * (b(i+1) - 2*b(i) + b(i-1)))... % TE_source
      + ((phi(i+2)-3*phi(i+1)+3*phi(i)-phi(i-1))/(24*dx(1)))... % TE_e
      - ((phi(i+1)-3*phi(i)+4*phi(i-1)-2*f0)/(24*dx(1))); % TE_w

i = n-1;
TE(i) = (dx(1)/24 * (b(i+1) - 2*b(i) + b(i-1)))... % TE_source
      + ((2*fn-4*phi(i+1)+3*phi(i)-phi(i-1))/(24*dx(1)))... % TE_e
      - ((phi(i+1)-3*phi(i)+3*phi(i-1)-phi(i-2))/(24*dx(1))); % TE_w

b= b_rb;
end;%diffusion1d




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
