function[DG] = DG_COMP(P,PL,PR,RHOL,RHOR,gam)

DG = G_DIFF(P,PL,gam,RHOL) + G_DIFF(P,PR,gam,RHOR);

end