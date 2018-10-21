function[PSTAR,USTAR] = STARVAL(PL,PR,UL,UR,RHOL,RHOR,gam)

P0 = 0.5*(PL + PR);
TOL = eps;
ERR = 1;
P_OLD = P0;
i = 0;
maxiter = 30;

while ERR > TOL && i<maxiter
    
    P = P_OLD - (G_COMP(P_OLD,PL,PR,UL,UR,RHOL,RHOR,gam)/DG_COMP(P_OLD,PL,PR,RHOL,RHOR,gam));
    ERR = abs(P - P_OLD);
    P_OLD = P;
    i = i + 1;
    
end

PSTAR = P;
USTAR = 0.5*(UL + UR) + 0.5*(-G_FUNC(PSTAR,PL,gam,RHOL) + G_FUNC(PSTAR,PR,gam,RHOR));

end