%% CEA Output Plotting for Mixture Ratio Selection
clear all; close all; clc;

% addpath("GCH4_LOX_PcSweep/")

names = [];
titles = [];
for i=1:15
    names = [names sprintf("GCH4_LOX_PcSweep/GCH4_LOX_PC%s",num2str(i+9))];
    titles = [titles sprintf("GCH4 & LO2, Pc = %s bar",num2str(i+9))];
end
filenames = strcat(names,".html");
imgnames = strcat(names,".png");
csvnames = strcat(names,".csv");

% Engine Parameters
Dc = 3.5*0.0254;% m
Ac = pi*Dc^2/4;
At = Ac/9;
Dt = sqrt(4/pi*At);

for fidx = 1:length(filenames)
    file = fopen(filenames(fidx));
    str = fscanf(file,'%s');
    
    % Find Chamber Pressure
    PCidx = strfind(str,"p,bar=")+6;
    Pc = str2num( str(PCidx:PCidx+strfind(str(PCidx:end),"#")-2) )*1e5;% Pa

    % Find O/F ratios used
    OFidx = strfind(str,"o/f=")+4;
    OF = str2num( str(OFidx:OFidx+strfind(str(OFidx:end),"#")-2) );
    n_OF = length(OF);

    % Find data
    EQidx = strfind(str,"PHI,EQ.RATIO=")+13;
    TCidx = strfind(str,"T,K")+3;
    CSidx = strfind(str,"CSTAR")+11;
    ISidx = strfind(str,"Isp")+9;
    ARidx = strfind(str,"Ae/At")+11;
    CFidx = strfind(str,"CF")+8;

    % Extract data
    arr = zeros(n_OF,11);
    arr(:,1) = OF;
    for i=1:n_OF
        arr(i,2) = str2num(str(EQidx(i):EQidx(i)+7));
        arr(i,3) = str2num(str(TCidx(i):TCidx(i)+6));
        arr(i,4) = str2num(str(CSidx(i):CSidx(i)+5));
        Ispidx = strfind(str(ISidx(i):ISidx(i)+5),".")+ISidx(i)+1;
        arr(i,5) = str2num(str(Ispidx:Ispidx+5));
        arr(i,6) = str2num(str(ARidx(i):ARidx(i)+5));
        arr(i,7) = str2num(str(CFidx(i):CFidx(i)+5));
        [arr(i,8), arr(i,9), arr(i,10), arr(i,11)] = get_thrust_mdot(At,Pc,OF(i),arr(i,7),arr(i,4));
    end

    % Convert data array to table
    col_labels = {'OF' 'EqRatio' 'Tc' 'Cstar' 'Isp' 'AeAt' 'CF' 'T' 'mdot' 'mdotO' 'mdotF'};
    tab = array2table(arr,'VariableNames',col_labels);
    tab.Properties.VariableUnits = {'' '' 'K' 'm/s' 'm/s' '' '' 'N' 'kg/s' 'kg/s' 'kg/s'};

    % Find maxima of data
    [Tc_m,Tc_idx] = max(tab.Tc);
    [Cstar_m,Cstar_idx] = max(tab.Cstar);
    [Isp_m,Isp_idx] = max(tab.Isp);

    % Plot data vs O/F
    figure
    yyaxis right
    plot(tab.OF,tab.EqRatio,'DisplayName','Eq. Ratio','LineWidth',2)
    hold on
    plot(tab.OF,tab.AeAt,'DisplayName','Ae/At','LineWidth',2)
    yyaxis left
    plot(tab.OF,tab.Tc,'DisplayName','Tc, K','LineWidth',2)
    plot(tab.OF,tab.Cstar,'DisplayName','C*, m/s','LineWidth',2)
    plot(tab.OF,tab.Isp,'DisplayName','Isp, m/s','LineWidth',2)
    scatter(OF(Tc_idx),Tc_m,'DisplayName',['Tc=' num2str(Tc_m) ', O/F=' num2str(OF(Tc_idx))])
    scatter(OF(Cstar_idx),Cstar_m,'DisplayName',['C*=' num2str(Cstar_m) ', O/F=' num2str(OF(Cstar_idx))])
    scatter(OF(Isp_idx),Isp_m,'DisplayName',['Isp=' num2str(Isp_m) ', O/F=' num2str(OF(Isp_idx))])
    legend('Location', 'east')
    xlabel('O/F Ratio')
    title(titles(fidx))
    hold off
    saveas(gcf, imgnames(fidx))
    close(gcf)
    
    % Print Table
%     csv = fopen(csvnames(fidx),"w");
%     fprintf('\n\n<strong>%s</strong>\nChamber ID: %0.4f in\nThroat ID: %0.4f in\nAc/At: %0.1f\nChamber Pressure: %0.1f bar\n\n',titles(fidx),Dc/0.0254,Dt/0.0254,Ac/At,Pc/1e5)
    writetable(tab,csvnames(fidx))
    %     disp(tab)
end

function [T, mdot, mdotO, mdotF] = get_thrust_mdot(At,Pc,OF,Cf,Cstar)
    T = Cf*Pc*At;
    mdot = Pc*At/Cstar;
    mdotO = mdot/(1 + 1/OF);
    mdotF = mdot - mdotO;
end

% function out = getSLPM(tab,OF)
% idx = find(tab.OF == OF);
% syms x
% out = vpasolve(SLPM2Mdot(16.043,x) == tab.mdotF(idx), x);
% end