function S = SatO2_2(C,Beta_p);
P = C/Beta_p;
S = (1 + 23400./((P).^3 + 150*P)).^(-1);
end