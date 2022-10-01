function CollectData
% % passiveStressInGallbladder(PatientNumber, num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones, ...
% %                                                                             ... rho, s1, s2, s3, plt)
% % Dat = [PatientNumber, PercentageOfLesion, D1, D2, D3, hGB, sigma_max_p];
% 
num = xlsread('CombinedData.xlsx');  % Read Table 1 of paper
% num  % needed for QA in num

Tot1 = [];
Tot2 = [];
dada = [];
rho = 0; % density of gallstones
Cases = 51;

for PercentChangeInhGB = -50:25:50
    for K = 1:Cases
%                           PatientNumber, num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones   rho, s1, s2, s3, plt
        Dat1 = passiveStressInGallbladder(K, num, PercentChangeInhGB,   66,                  0,                0,  2, 2, 5,     0); % 1/3 part, (66% cut out), with higher hGB
        Tot1 = [Tot1; Dat1(7)];
        hGB_contracted = Dat1(6);
        Dat2 = passiveStressInGallbladder(K, num, 0,                    33,                  0,                0,  1, 1, 1,      0); %2/3 part, (33% cut out), with same hGB
        Tot2 = [Tot2; Dat2(7)];
    end
    dada = [dada; hGB_contracted mean(Tot1) std(Tot1) min(Tot1) max(Tot1) mean(Tot2) std(Tot2) min(Tot2) max(Tot2)];
%    [size(Tot1) size(Tot2) size(dada)]
    Tot1 = [];
    Tot2 = [];              
end 

xlswrite('ModelCompare.xls', dada)

 
