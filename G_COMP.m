function[G] = G_COMP(P,PL,PR,UL,UR,RHOL,RHOR,gam)

DU = UR - UL;
[GL,~,~,~] = G_FUNC(P,PL,gam,RHOL);
[GR,~,~,~] = G_FUNC(P,PR,gam,RHOR);
G = GL + GR + DU;

end