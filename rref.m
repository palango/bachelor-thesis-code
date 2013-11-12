function [X2] = rref(N,X,ERR)

%N = 10;
%X = linspace(0,1,N+1);
XC = (X(1:N)+X(2:N+1))/2;
DX = X(2:N+1)-X(1:N);
%ERR = cos(pi*XC);
ERR=abs(ERR);

SUMERR=0;
for I=1:N
  SUMERR = SUMERR + DX(I)*ERR(I);
end
EERR = SUMERR / N;

X2 = zeros(1,N+1);
REST=0;
XCUR = 0;
IDX=2;
for I=1:N
  XBCELL = X(I);
  XCUR = XBCELL;
  while 1
    % how much's missing?
    MISSING = EERR - REST;
    RESTC = ERR(I) * (DX(I) - (XCUR-XBCELL));
    % can call be filled?
    if RESTC >= MISSING
      DELTAX = MISSING/ERR(I);
      XCUR=XCUR+DELTAX;
      X2(IDX) = XCUR;

      % get ready for next iteration
      MISSING=EERR;
      REST=0;
      IDX=IDX+1;
    % take rest of cell an jump to next
    else
      REST=REST+RESTC;
      break;
    end
  end
end

end %function
