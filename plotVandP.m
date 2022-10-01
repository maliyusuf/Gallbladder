function plotVandP

GallBladderData  % Loads the GallBladderData.m file

gamm = Ve - C * (pe - pd);
R = (texp - ti) / (C * log((V0-gamm)/(Vexp - gamm)));  % Flow resistance
te = log((V0-gamm)/(Ve-gamm)) / log((V0-gamm)/(Vexp-gamm)) * (texp - ti);
t = ti:(te+ti);
press = pd +   (pe - pd) * exp( (te+ti-t) / (R*C)  );
Vol  = C  * ( (pe - pd) * exp( (te+ti-t) / (R*C)  ) ) + gamm;

% Start Plotting

scrsz = get(0,'ScreenSize'); % [left, bottom, width, height]:
figure('OuterPosition',[1 5 scrsz(3) scrsz(4)])

subplot(2,1,1)
    plot(t,Vol, 'LineWidth',2)
    fsize = 16;
    h1=xlabel('Time (min)');
    h2=ylabel('Volume (mL)');
    title(['te = ', num2str(te)],'FontSize',fsize );
    set(gca,'fontsize',fsize) % increase the size
    set(h1,'fontsize',fsize) % increase the size
    set(h2,'fontsize',fsize) % increase the size	
    grid on

subplot(2,1,2)
    plot(t,press, 'LineWidth',2)
    fsize = 16;
    h1=xlabel('Time (min)');
    h2=ylabel('Pressure (mmHg)');
    set(gca,'fontsize',fsize) % increase the size
    set(h1,'fontsize',fsize) % increase the size
    set(h2,'fontsize',fsize) % increase the size	
    grid on

