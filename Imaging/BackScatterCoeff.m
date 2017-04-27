clc ;
clear all;

% calculation of back scatter coefficient (Love & Scott , 1978)
% created 5/14/2013 by Ashwin Vishwanathan

% assuming weight factor of Os>Ua>Pb

WOs = 0.6;
Wua = 0.3;
Wpb = 0.1;

ZOs = 76;
Zua = 92;
Zpb = 82;

%E = 10 ; % in KeV
E = 1:0.25:10;

% for i = 1:1:length(WOs)
% 
% Zavg(i) = WOs(i)*ZOs + Wua*Zua + Wpb*Zpb;
% nbtwenty(i) = -5.23791e-3 + (1.5048371e-2*Zavg(i)) - (1.67373e-4*Zavg(i)^2) + (7.16e-7*Zavg(i)^3);
% a(i) = -0.11128 + (3.0289e-3*Zavg(i)) - (1.5498e-5*Zavg(i)^2) ;
% 
% end

Zavg = WOs*ZOs + Wua*Zua + Wpb*Zpb;

nbtwenty = -5.23791e-3 + (1.5048371e-2*Zavg) - (1.67373e-4*Zavg^2) + (7.16e-7*Zavg^3);

a = -0.11128 + (3.0289e-3*Zavg) - (1.5498e-5*Zavg^2) ;

clear i;



for i = 1:1:length(E)

n(i) = nbtwenty*(1+ a*log(E(i)/20)) ; 

end

%% 

plot (E,n)
xlabel ('E in Kev');
Ylabel ('Backscatter coefficient \eta ');



