#!/home/kmtu/bin/octave -qf

clear all;

if (nargin() < 1)
    error("Usage: $plotSD.m <sdresult>")
endif

sd_file = argv(){1}
%set(0, "defaultaxesfontname", "Arial")
%set(0, "defaulttextfontname", "Arial") 
set(0, "defaultlinelinewidth", 4);
sd = load(sd_file);

%maxlag = size(sd.timeLags, 1)
maxlag = 1000

f = figure(1);
clf;
%set(f, "paperorientation", "portrait")

subplot(2, 1, 1)
[c,h_contourf] = contourf(sd.rBins, sd.timeLags(1:maxlag), sd.sdCorr(1:maxlag,:,1),5);
title(strcat("Correlation 1 - lag 1:", num2str(maxlag)));
xlabel("r (nm)");
ylabel("lag (ps)");

subplot(2, 1, 2)
h_mesh = mesh(sd.rBins, sd.timeLags(1:maxlag), sd.sdCorr(1:maxlag,:,1));
xlabel("r (nm)");
ylabel("lag (ps)");

%print(strcat('sdCorr1', '-E4E4-S1-', num2str(maxlag), '.eps'), '-deps', '-color', '-S800,1200');
print(strcat('sdCorr1', '-E6E5-S2-', num2str(maxlag), '.png'), '-dpng', '-color', '-S800,1200');


%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

%legend("{/Symbol D}t = 1 fs");
%legend("Total", "Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Na", "Cross Cl-Cl", "location", "northwest");
%axis([0,100,-600, 1000]);


%axis([0,0.1,0,1.5]);
%print('ecNoAverageCesaro-zoom.eps', '-deps', '-color');

