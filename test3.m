num = xlsread('CombinedData.xlsx');  % Read Table 1 of paper

Tot1 = [];
Tot2 = [];
dada = [];

PercentChangeInhGB = 0;

for K = 1:88
     Dat1 = passiveStressInGallbladder(K, num,PercentChangeInhGB,0, 0, 0, 1,1,1, 0); % 1/3 part, (66% cut out), with higher hGB
     Tot1 = [Tot1; Dat1(7)];
end
Tot1 