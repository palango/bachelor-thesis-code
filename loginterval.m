% Vector von logarithmisch verteilen Elementen in [from, to] erzeugen
function x = loginterval(from, to, n)
x = from + (logspace(0, log10(11), n)-1)./(10/(to-from));
end;
