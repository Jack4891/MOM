function GFunc=GFuncNonSingPart(R,k0)

GFunc=(exp(-1i.*k0.*R)./R)-(1./R);
GFunc(isnan(GFunc))=-1i*k0;