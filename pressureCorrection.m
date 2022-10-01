function pressureCorrection
% function pressureCorrection
%
%  PercentOfGallstones  = percentage of the total gallbladder volume
%  rho                  = density of gallstones in Kg/m^3
%  s1, s2, s3           = scaling factors for the contracted gallbladder 

num = xlsread('data2010.xlsx');  % Read Table 1 of paper by Li 2010

granddada = [];

sa = [1 1 1];   % Dimensions scaling factors
sb = [1/2 1/2 1/5];   % Dimensions scaling factors

dada = [];
for PercentOfGallstones = 33:33:99
    dat = go(num, PercentOfGallstones, sa);
    dada = [dada; dat];
end
granddada = [granddada; dada];

dada = [];
for PercentOfGallstones = 33:33:99
    dat = go(num, PercentOfGallstones, sb);
    dada = [dada; dat];
end
granddada = [granddada; dada];

  xlswrite('PressureCorrection.xls', granddada)
end 

function  dat = go(num, PercentOfGallstones, s)
        dat = [];

            for rho = 1000:1000:5000
                datum = generate(num, PercentOfGallstones, rho, s);
                dat = [dat; datum];
            end
end 
function datum = generate(num, PercentOfGallstones, rho, s)
        report = [];
        datum = [];

            for I = 1:51
               p_correction = call_correction(I, num, s, PercentOfGallstones, rho);
               report = [report; p_correction];
            end
            
            [C1 L1] = min(report);
            [C2 L2] = max(report);
        datum = [datum; PercentOfGallstones, rho, s , L1, C1, L2, C2];
       
    end 
function p_correction = call_correction(I, num, s, PercentOfGallstones, rho)

pe = 11;    % gallbladder internal pressure while emptying, 11 mmHg in the paper
            % Ref: Majeed - Continuous in vivo manometry of the human gallbladder - 1994

D3 = num(I,4);  % Diameter along the major axis at time t, in mm
k1 = num(I,2);  % k1 = D1/D3 ;    
k2 = num(I,3);  % k2 = D2/D3;
D1 = k1*D3;
D2 = k2*D3;

g = 9.8;       % m/sec^2

a = D1/2/1000;        % Convert mm to meters tby dividing with 1000
b = D2/2/1000;
c = D3/2/1000;

MinimumSurfaceArea = pi*a*b;            % This is what you see from the top, even if there is a contraction 
Volume_of_ellipsoid = (4/3)*a*b*c;      % Volume BEFORE contraction 
frac = PercentOfGallstones/100;         % Fraction of gallstone volume in the gallbladder
VolumeOfStones = frac*Volume_of_ellipsoid;        % Volume of gallstones in mm^3
F = rho*g*VolumeOfStones;               % Lengths measured in mm, Force in Newtons 
p = F / MinimumSurfaceArea / 133 ; % Divide by two as stones area on one side only
                                            %133 Converts pascals into mmHg
p_correction = p *100 / pe; 
                                            
end