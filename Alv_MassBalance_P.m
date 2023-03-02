function P = Alv_MassBalance_P(V2,V1,t2,t1,PcO2,PAO21,pars)
% Physiological parameters
    yO2 = pars(1); 
    PB = pars(2); 
    PH2O = pars(3); 
    R = pars(4);
    T = pars(5);
    Beta_p = pars(6);
    k = pars(7);
     
% Alveolar mass balance
PIO2 = yO2*(PB-PH2O);
dt = t2-t1;
fact1 = (V2-V1)./V1;
fact2 = (V2-V1)/dt;
coeff = k*R*T*dt*Beta_p./V1;
% Inspiration
num1 = fact1.*(PIO2-PAO21) - coeff.*(PAO21-PcO2) + PAO21;

%Expiration
num2 = -coeff.*(PAO21-PcO2) + PAO21;
    
    %Inspiration
    if (fact2 >= 0)
        P = num1;
    else
    %Expiration
        P =  num2;
    end
end