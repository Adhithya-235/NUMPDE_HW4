function[Gpd] = G_DIFF(P,Pd,gam,rhod)

Cd = sqrt(gam*Pd/rhod);
Ad = 2/((gam + 1)*rhod);
Bd = Pd*((gam - 1)/(gam + 1));


if P>Pd
    
    Gpd = sqrt(Ad/(Bd + P))*(1 - ((P - Pd)/(2*(Bd + P))));
    
elseif P<=Pd
    
    Gpd = (1/(rhod*Cd))*(((P/Pd)^((-1 - gam))/(2*gam)));
    
end


end