function stressAndLesion

 PUpper = [];
 sUpper = [];
 PLower = [];
 sLower = [];
num = xlsread('CombinedData.xlsx');  % Read Table 1 of paper

for PercentageOfLesion = 10:5:90
% % passiveStressInGallbladder(PatientNumber, num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones, ...
% %                                                                             ... s1, s2, s3, plt)
% % Dat = [PatientNumber, PercentageOfLesion, D1, D2, D3, hGB, sigma_max_p];
 DatUpper = passiveStressInGallbladder(9, num, 0, PercentageOfLesion, 0, 1, 1,1, 0);
 PUpper = [PUpper; PercentageOfLesion];
 sUpper = [sUpper; DatUpper(7)];     % DatUpper(7) has the value of hGB recorded
 DatLower = passiveStressInGallbladder(9, num, 0, 100-PercentageOfLesion, 0,1, 1,1, 0);
 PLower = [PLower; PercentageOfLesion];
 sLower = [sLower; DatLower(7)];
 
end

fsize = 16;
plot(PUpper,sUpper,'b*-')
hold on
plot(PLower,sLower,'rx-', 'LineWidth',2)
h1 = xlabel(['Percentage of Lesion']);
h2 = ylabel(['max(\sigma) mmHg']);
title(['Upper h_G_B = ',num2str(DatUpper(6)), ' mm', ', Lower h_G_B = ',num2str(DatLower(6)),' mm'],'FontSize',fsize)
set(gca,'fontsize',fsize) % increase the size
set(h1,'fontsize',fsize) % increase the size
set(h2,'fontsize',fsize) % increase the size	
legend('Upper Part','Lower Part')

pause;
set(0,'ShowHiddenHandles','on')
delete(get(0,'Children'))
        
 