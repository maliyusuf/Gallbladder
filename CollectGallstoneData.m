function CollectGallstoneData
% % passiveStressInGallbladder(PatientNumber, num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones,... 
%                                                                                 ... rho,s1, s2, s3, plt)
% % Dat = [PatientNumber, PercentageOfLesion, D1, D2, D3, hGB, sigma_max_p];
% 
num = xlsread('CombinedData.xlsx');  % Read Table 1 of paper
% num  % needed for QA in num

Tot2 = [];
dada = [];

for rho = 1000:1000:5000
    for K = 1:51       
        PercentChangeInhGB = 0;
        PercentageOfLesion = 0;
        PercentOfGallstones = 33;
        Dat2 = passiveStressInGallbladder(K,num,PercentChangeInhGB,PercentageOfLesion, PercentOfGallstones, rho, 1, 1, 1, 0); 		%2/3 part, (33% cut out), with same hGB
        Tot2 = [Tot2; Dat2(7)];
    end
    dada = [dada; 33 rho mean(Tot2) std(Tot2) min(Tot2) max(Tot2)];
    Tot2 = [];
end

Tot2 = [];
    
for rho = 1000:1000:5000
    for K = 1:51
        PercentChangeInhGB = 0;
        PercentageOfLesion = 0;
        PercentOfGallstones = 66;
        Dat2 = passiveStressInGallbladder(K,num,PercentChangeInhGB,PercentageOfLesion, PercentOfGallstones, rho, 1, 1, 1, 0); 		%2/3 part, (33% cut out), with same hGB
        Tot2 = [Tot2; Dat2(7)];
    end
    dada = [dada; 66 rho mean(Tot2) std(Tot2) min(Tot2) max(Tot2)];
    Tot2 = [];
end

xlswrite('ModelwithGallstones.xls', dada)

end 

