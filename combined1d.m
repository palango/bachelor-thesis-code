function [xc, phi, dx, A, b, TE] = combined1d(n, f0, fn, f, v, rho, alpha, beta)
p=1;q=1;
x = linspace(0, 1, n + 1);
xc = zeros(n, 1);
xc(1:n) = 0.5 * (x(1:n) + x(2:n+1));

dx = x(2:n+1) - x(1:n);

% Gleichungssystem für äquidistantes Gitter
A = zeros(n);
a = rho * v / (dx(1));
c = alpha / dx(1)^2;

A = diag(vse(q*a*(1-beta)-p*2*c, n)) + diag(vse(q*a*beta/2+p*c, n-1), 1) + diag(vse(q*(-a)*(1-beta/2)+p*c, n-1), -1);

% Quellterme in Mittelpunkten berechnen
b = arrayfun(f, xc);

% Randbedingungen in A
A(1, 1) = q*a*(1-beta/2)-3*c*p;
A(n, n) = q*(-a)*(beta/2)-p*3*c;
% Randbedingungen in b einfügen
b(1) = b(1) + a*f0*q-p*2*c*f0;
b(n) = b(n) - a*fn*q-p*2*c*fn;

phi = pinv(A)*b;

% Truncation Error berechnen
for i=3:n-2
  TE(i) = (dx(1)/24 * (b(i+1) - 2*b(i) + b(i-1)))... % TE_source
        + ((phi(i+2)-3*phi(i+1)+3*phi(i)-phi(i-1))/(24*dx(1)))... % TE_e
        - ((phi(i+1)-3*phi(i)+3*phi(i-1)-phi(i-2))/(24*dx(1))); % TE_w
end;
end
