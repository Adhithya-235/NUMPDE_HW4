function euler_solver
    PL = 1;
    PR = 0.1;
    UL = 0;
    UR = 0;
    RHOL = 1;
    RHOR = 0.125;
    gam = 1.4;
    T = .25;
    [~,CL,~,~] = G_FUNC(0,PL,gam,RHOL);
    [~,CR,~,~] = G_FUNC(0,PR,gam,RHOR);
    N = 1000;
    x = linspace(-0.5,0.5,N);
    S = x/T;
    [PSTAR,USTAR] = STARVAL(PL,PR,UL,UR,RHOL,RHOR,gam);
    rho_vec = zeros(1,N);
    u_vec = rho_vec;
    p_vec = rho_vec;
    
    
    for i=1:length(S)
       
        if(S(i) < USTAR)
           %we are to the left of the contact wave
           if(PSTAR > PL)
              %we have a shock to the left of the contact wave
              type = 'shock';
              [SL,~,~,RHOSTAR] = speed_calc_left(UL,PSTAR,PL,gam,...
                                     USTAR,CL,RHOL,type);
              if(S(i) < SL)
                 %we are left of the shock
                 u_vec(i) = UL;
                 rho_vec(i) = RHOL;
                 p_vec(i) = PL;
              else
                 u_vec(i) = USTAR;
                 rho_vec(i) = RHOSTAR;
                 p_vec(i) = PSTAR;
              end
           else
              %we are in a fan
              type = 'rarefaction';
              [~,SHL,STL,RHOSTAR] = speed_calc_left(UL,PSTAR,PL,gam,...
                                     USTAR,CL,RHOL,type);
              if(S(i) < SHL)
                 %we are left of the fan
                 u_vec(i) = UL;
                 rho_vec(i) = RHOL;
                 p_vec(i) = PL;
              elseif(S(i) < STL)
                 %we are in the fan
                 [RHOF,UF,PF] = wfan(RHOL,gam,CL,S(i),PL,UL);
                 u_vec(i) = UF;
                 rho_vec(i) = RHOF;
                 p_vec(i) = PF;
              else
                 % we are right of the fan
                 u_vec(i) = USTAR;
                 rho_vec(i) = RHOSTAR;
                 p_vec(i) = PSTAR;  
              end
           end
        elseif(S(i) == USTAR)
            %do nothing
        else
           %we are to the right of the contact wave
           if(PSTAR > PR)
              %we have a shock to the right of the contact wave
              type = 'shock';
              [SR,~,~,RHOSTAR] = speed_calc_right(UR,PSTAR,PR,gam,...
                                     USTAR,CR,RHOR,type);
              if(S(i) > SR)
                 %we are right of the shock
                 u_vec(i) = UR;
                 rho_vec(i) = RHOR;
                 p_vec(i) = PR;
              else
                 u_vec(i) = USTAR;
                 rho_vec(i) = RHOSTAR;
                 p_vec(i) = PSTAR;
              end
           else
              %we are in a right fan
              type = 'rarefaction';
              [~,SHR,STR,RHOSTAR] = speed_calc_right(UL,PSTAR,PR,gam,...
                                     USTAR,CR,RHOR,type);
              if(S(i) < SHR)
                 %we are right of the fan
                 u_vec(i) = UR;
                 rho_vec(i) = RHOR;
                 p_vec(i) = PR;
              elseif(S(i) < STR)
                 %we are in the fan
                 [RHOF,UF,PF] = wfan(RHOR,gam,-CR,S(i),PR,UR);
                 u_vec(i) = UF;
                 rho_vec(i) = RHOF;
                 p_vec(i) = PF;
              else
                 % we are left of the fan
                 u_vec(i) = USTAR;
                 rho_vec(i) = RHOSTAR;
                 p_vec(i) = PSTAR;  
              end
           end   
        end
    end
    figure(1)
    title('Primitive Variables at t = 0.25');
    subplot(2,2,1)
    plot(x,rho_vec)
    ylabel('$\rho$','interpreter','latex')
    xlabel('$x$','interpreter','latex')
    subplot(2,2,2)
    plot(x,u_vec)
    ylabel('$u$','interpreter','latex')
    xlabel('$x$','interpreter','latex')
    subplot(2,2,3)
    plot(x,p_vec)
    ylabel('$p$','interpreter','latex')
    xlabel('$x$','interpreter','latex')
    subplot(2,2,4)
    E = p_vec/(gam-1) + 0.5*rho_vec.*u_vec.^2;
    plot(x,E)
    ylabel('$E$','interpreter','latex')
    xlabel('$x$','interpreter','latex')
return

function [SL,SHL,STL,RHOSTAR]=speed_calc_left(UL,PSTAR,PL,gam,USTAR,CL,...
                              RHOL,type)
    if(strcmp(type,'shock'))
        val = (gam+1)*PSTAR/(2*gam*PL);
        SL = UL - CL*sqrt(val + (gam-1)/(2*gam));
        SHL = -1;
        STL = -1;
        denom = (gam-1)/(gam+1)*(PSTAR/PL) + 1;
        RHOSTAR = RHOL*((gam-1)/(gam+1) + PSTAR/PL)/denom;
    else
        %rarefaction wave
        SL = -1;
        SHL = UL - CL;
        STL = USTAR - CL*power(PSTAR/PL,(gam-1)/(2*gam));
        RHOSTAR = RHOL*power(PSTAR/PL,1/gam);
    end
return

function [SR,SHR,STR,RHOSTAR]=speed_calc_right(UR,PSTAR,PR,gam,USTAR,CR,...
                              RHOR,type)
    if(strcmp(type,'shock'))
        val = (gam+1)*PSTAR/(2*gam*PR);
        SR = UR + CR*sqrt(val + (gam-1)/(2*gam));
        SHR = -1;
        STR = -1;
        denom = (gam-1)/(gam+1)*(PSTAR/PR) + 1;
        RHOSTAR = RHOR*((gam-1)/(gam+1) + PSTAR/PR)/denom;
    else
        %rarefaction wave
        SR = -1;
        SHR = UR + CR;
        STR = USTAR + CR*power(PSTAR/PR,(gam-1)/(2*gam)); 
        RHOSTAR = RHOR*power(PSTAR/PR,1/gam);
    end
return

function [RHOF,UF,PF] = wfan(rho,gam,C,S,P,U)
    %NOTE: If we are using the right expression, must input
    %-CR for C!!!!!
    val = 2/(gam+1) + ((gam-1)/((gam+1)*C))*(U-S);
    RHOF = rho*power(val,2/(gam-1));
    UF = 2*(C + 0.5*(gam-1)*U + S)/(gam+1);
    PF = P*power(val,2*gam/(gam-1));
return