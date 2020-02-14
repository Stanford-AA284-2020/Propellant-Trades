%% Chamber Pressure Trends

% addpath("GCH4_LOX_PcSweep/")

names = [];
% titles = [];
for i=1:15
    names = [names sprintf("GCH4_LOX_PcSweep/GCH4_LOX_PC%s",num2str(i+9))];
%     titles = [titles sprintf("GCH4 & LO2, Pc = %s bar",num2str(i+9))];
end
% filenames = strcat(names,".html");
% imgnames = strcat(names,".png");
csvnames = strcat(names,".csv");

storetab = readtable(csvnames(1));
storetab = storetab(1,:);
storetab{1,:} = NaN;
Pc = zeros(length(names),1);
for i=1:15
    Pc(i) = str2double(names{i}(end-1:end));
    tab = readtable(csvnames(i));
    storetab = [storetab;tab(11,:)];
end
storetab = rmmissing(storetab);
storetab = addvars(storetab,Pc,'Before','OF');
plot(storetab.Pc,storetab.Cstar)