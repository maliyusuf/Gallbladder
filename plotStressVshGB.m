function plotStressVshGB

% Dat = passiveStressInGallbladder(PatientNumber, num, PercentChangeInhGB, PercentageOfLesion, PercentOfGallstones, rho, s1, s2, s3, plt)
% Dat = [PatientNumber, PercentageOfLesion, D1, D2, D3, hGB, sigma_max_p];

num = xlsread('CombinedData.xlsx');  % Read Table 1 of paper
A = [];
B = [];
Cases = 51;
num = num(1:51,:);
rho = 0;   % density of gallstones 

PercentChangeInhGB = -50;
    for I = 1:Cases
        Dat = passiveStressInGallbladder(I, num, PercentChangeInhGB, 50, 0, rho, 1, 1, 1, 0); 
        A = [A; Dat(6) Dat(7)];           % Dat6 = hGB, Dat7 = stress
    end
    B = [B; Dat(6) mean(A(:,2)) std(A(:,2)) min(A(:,2)) max(A(:,2))];
    A = [];

PercentChangeInhGB = -25;
    for I = 1:Cases
        Dat = passiveStressInGallbladder(I, num, PercentChangeInhGB, 0, 0, rho, 1, 1, 1, 0); 
        A = [A; Dat(6) Dat(7)];
    end
    B = [B; Dat(6) mean(A(:,2)) std(A(:,2)) min(A(:,2)) max(A(:,2))];
    A = [];
    
PercentChangeInhGB = 0;
    for I = 1:Cases
        Dat = passiveStressInGallbladder(I, num, PercentChangeInhGB, 0, 0, rho, 1, 1, 1, 0); 
        A = [A; Dat(6) Dat(7)];
    end
    B = [B; Dat(6) mean(A(:,2)) std(A(:,2)) min(A(:,2)) max(A(:,2))];
    A = [];
    
PercentChangeInhGB = 25;
    for I = 1:Cases
        Dat = passiveStressInGallbladder(I, num, PercentChangeInhGB, 0, 0, rho, 1, 1, 1, 0); 
        A = [A; Dat(6) Dat(7)];
    end
    B = [B; Dat(6) mean(A(:,2)) std(A(:,2)) min(A(:,2)) max(A(:,2))];
    A = [];
        
PercentChangeInhGB = 50;
    for I = 1:Cases
        Dat = passiveStressInGallbladder(I, num, PercentChangeInhGB, 0, 0, rho, 1, 1, 1, 0); 
        A = [A; Dat(6) Dat(7)];
    end
    B = [B; Dat(6) mean(A(:,2)) std(A(:,2)) min(A(:,2)) max(A(:,2))];
    A = [];

xlswrite('StressVshGB.xls', B)

% scrsz = get(0,'ScreenSize');
% figure('OuterPosition',[1 5 scrsz(3) scrsz(4)]);  
% plot(A(:,1),A(:,2),'r-','LineWidth',2)
% h1 = xlabel('h_G_B (mm)');
% h2 = ylabel('Max(\sigma) (mmHg)');
% fsize = 20;
% set(gca,'fontsize',fsize) % increase the size
% set(h1,'fontsize',fsize) % increase the size
% set(h2,'fontsize',fsize) % increase the size

% pause;
% set(0,'ShowHiddenHandles','on')
% delete(get(0,'Children'))

end 
