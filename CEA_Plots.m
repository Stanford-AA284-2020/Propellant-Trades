%% CEA Output Plotting for Mixture Ratio Selection
clear all; close all; clc;

% filenames = ["CEA_CH4_OX_MRrange144.html","CEA_LNG_OX_MRrange.html",...
%     "CEA_LNGworst_OX_MRrange.html"];
% filenames = ["CEA_CH4_OX_MRrange_Shifting.html","CEA_LNG_OX_MRrange_Shifting.html",...
%     "CEA_LNGworst_OX_MRrange_Shifting.html"];
% titles = ["CH4(L) & O2(L)",...
%     "LNG [92% CH4(L), 7%C2H6(L), 1%C3H8(L)] & O2(L)",...
%     "LNG [86% CH4(L), 11%C2H6(L), 3%C3H8(L)] & O2(L)"];
% filenames = ["CEA_RP1_OX_MRrange144.html"];
% titles = ["RP-1 & O2(L)"];
% filenames = ["CEA_RP1_N2O_MRrange.html"];
% titles = ["RP-1 & N2O"];
% filenames = ["CEA_LNGworst_OX_MRrange.html","CEA_RP1_OX_MRrange144.html","CEA_C3H8_OX_MRrange.html"];
% titles = ["LNG [86% CH4(L), 11%C2H6(L), 3%C3H8(L)] & O2(L)","RP-1 & O2(L)","C3H8(L) & O2(L)"];
% filenames = ["CEA_C3H8_N2O_MRrange.html"];
% titles = ["C3H8(L) & N2O"];

% filenames = ["CEA_CH4_OX_MRrange144.html"];
% titles = ["CH4(L) & O2(L)"];
% filenames = ["CEA_GCH4_O2_MRrange.html"];
% titles = ["CH4 & O2(L), Pc = 24 bar"];
% filenames = ["CEA_CH4_OX_14bar.html"];
% titles = ["CH4 & O2, Pc = 14 bar"];
filenames = ["CEA_CH4_OX_10bar.html"];
titles = ["CH4 & O2, Pc = 10 bar"];

% Engine Parameters
Dc = 3.5*0.0254;% m
Ac = pi*Dc^2/4;
At = Ac/9;
Dt = sqrt(4/pi*At);
% Pc = 1.4e6;% Pa

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
    
    % Print Table
    fprintf('\n\n<strong>%s</strong>\nChamber ID: %0.4f in\nThroat ID: %0.4f in\nAc/At: %0.1f\nChamber Pressure: %0.1f bar\n\n',titles(fidx),Dc/0.0254,Dt/0.0254,Ac/At,Pc/1e5)
    disp(tab)
%     FSLPM = getSLPM(tab, 2.6)
end

function [T, mdot, mdotO, mdotF] = get_thrust_mdot(At,Pc,OF,Cf,Cstar)
    T = Cf*Pc*At;
    mdot = Pc*At/Cstar;
    mdotO = mdot/(1 + 1/OF);
    mdotF = mdot - mdotO;
end

function out = getSLPM(tab,OF)
idx = find(tab.OF == OF);
syms x
out = vpasolve(SLPM2Mdot(16.043,x) == tab.mdotF(idx), x);
end