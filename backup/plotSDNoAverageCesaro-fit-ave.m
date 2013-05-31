#!/home/kmtu/bin/octave -qf

clear all;
global numIonTypes;

if (nargin() < 3)
    error("Usage: $plotNoAverageCesaro.m <dataBaseName> <skip> <dt>")
endif

dataBaseName = argv(){1}
skip = str2num(argv(){2}) #skipped interval in sdCorr data
deltaStep = str2num(argv(){3})

set(0, "defaultlinelinewidth", 4);

%md.dc = load("md-lag20000.dcNoAverageCesaro" );
%md3.dc = load("md3-lag20000.dcNoAverageCesaro" );

numMD = 1;

%initialize = load(strcat("./md0/", dataBaseName, num2str(deltaStep)));
%timeLags = initialize.timeLags;
%clear initialize;

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
    index = (idx1 - 1) * numIonTypes + idx2;
endfunction

# md(sdCorr_timeLag, sdCorr_rBin, sdCorrIonTypePairIndex, fileIndex)
for n = [1:numMD]
    data = load(strcat("./md", num2str(n-1), "/", dataBaseName, num2str(deltaStep)));
    md_tot(n) = data.ecSDTotalNoAverageCesaro;
    for i = [1:numIonTypes**2]
        md(:, :, :, n) = data.ecSDCorrNoAverageCesaro;
    endfor
endfor
timeLags = data.timeLags;
clear(data);

md_ave = squeeze(mean(md, 4))';
md_std = squeeze(std(md, 0, 4))';
md_err = squeeze(std(md, 0, 4) ./ sqrt(numMD))';  # standard error 
numIonTypePairs = size(md_ave, 3)
num_rBins = size(md_ave, 2)

fitRange = [20, 40; 40, 60; 60, 80; 80, 100]; #ps
fitRange *= floor(1000 / skip / deltaStep); #fs (frame)

# calculate slope for each segment of md_ave
for i = [1:numIonTypePairs]
  for j = [1:num_rBins]
    for r = [1:size(fitRange, 1)]
        slope(r,j,i) = polyfit(timeLags(fitRange(r, 1):fitRange(r, 2)), md_ave(fitRange(r, 1):fitRange(r, 2), j, i), 1)(1);
    endfor
  endfor
endfor

# evaluate the uncertainty in the slope of the fitting line
# reference: Numerical Recipes Chapter 15.2 (p.656)
for i = [1:numIonTypePairs] 
  for j = [1:num_rBins]
    for r = [1:size(fitRange, 1)]
        rec_sig2 = 1 ./ (md_std(fitRange(r, 1):fitRange(r, 2), j, i) .^ 2);
        S(r, j, i) = sum(rec_sig2, 1);
        Sx(r, j, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* rec_sig2, 1); 
        Sxx(r, j, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)).^2 .* rec_sig2, 1); 
        Sy(r, j, i) = sum(md_ave(fitRange(r, 1):fitRange(r, 2), j, i) .* rec_sig2, 1);
        Syy(r, j, i) = sum(md_ave(fitRange(r, 1):fitRange(r, 2), j, i).^2 .* rec_sig2, 1);
        Sxy(r, j, i) = sum(timeLags(fitRange(r, 1):fitRange(r, 2)) .* md_ave(fitRange(r, 1):fitRange(r, 2), j, i) .* rec_sig2, 1); 
    endfor
  endfor
endfor
Delta = S .* Sxx - Sx .* Sx;

# output slope for double check
slope_b = (S .* Sxy - Sx .* Sy) ./ Delta
slopeSD = sqrt(S ./ Delta);
# slopeSD(fitRange, corrIndex)

#numPlots = 1 + numIonTypes + numIonTypes*numIonTypes;

# standard error for selected values
errFrameInterval = floor(5000 / skip / deltaStep);
errFrames = [1:errFrameInterval:length(timeLags)]';


f1 = figure(1);
clf;
hold on;
errOffset = floor(200 / skip / deltaStep);

for i = [1:numIonTypePairs]
    drawFrames = errFrames .+ (i-1)*errOffset;
    drawFrames = drawFrames(drawFrames <= length(timeLags));
    if (i == 6)
#        p(i).plot = plot(timeLags, md_ave(:,i), "-", "color", [0, 0.6, 0]);
        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
        set(p(i).errorbar(1), "color", [0, 0.6, 0]);
    elseif (i == 7)
#        p(i).plot = plot(timeLags, md_ave(:,i), "-", "color", [0, 0, 0]);
        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
        set(p(i).errorbar(1), "color", [0, 0, 0]);
    else
        plotFormat = strcat("-", num2str(i));
        errorbarFormat = strcat("~", num2str(i));
        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i), errorbarFormat);
#        p(i).plot = plot(timeLags, md_ave(:,i), plotFormat);
    endif
endfor

legend("Total", "Na-Na", "Na-Cl", "Cl-Na", "Cl-Cl", "location", "northwest");

left = 25;
h_space = 20;
top = 972;
v_space = 44;
for i = [1:size(md_ave, 2)]
    for r = [1:size(fitRange, 1)]
        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(slope(r, i)));
    endfor
endfor

%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

%title(strcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution - md", "-ave"));
title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",\
             "data interval = ", num2str(skip), ", integration interval = ", num2str(deltaStep)));
%legend("{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Non-averaged partial sum (ps*S/m)");
axis([0,100,-600, 1000]);

print(strcat('ecNoAverageCesaro-', 'ave-skip-', num2str(skip), '-dt-', num2str(deltaStep), '.eps'), '-deps', '-color');
hold off
%axis([0,0.1,0,1.5]);
%print('ecNoAverageCesaro-zoom.eps', '-deps', '-color');

##############################

# standard deviation for selected values
for i = [1:size(fitRange, 1)]
    stdFrames(i) = (fitRange(i, 1) + fitRange(i, 2)) / 2;
endfor

f2 = figure(2);
%clf;
hold on;
stdOffset = floor(100 / skip / deltaStep);

for i = [1:numIonTypePairs]
    drawFrames = stdFrames .+ (i-1)*stdOffset;
    if (i == 6)
        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i));
        set(p(i).errorbar(1), "color", [0, 0.6, 0]);
    elseif (i == 7)
        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i));
        set(p(i).errorbar(1), "color", [0, 0, 0]);
    else
        plotFormat = strcat("-", num2str(i));
        errorbarFormat = strcat("~", num2str(i));
        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i), errorbarFormat);
    endif
endfor

legend("Total", "Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Na", "Cross Cl-Cl", "location", "northwest");

%left = 25;
%h_space = 20;
%top = 9.75;
%v_space = 0.41;
%for i = [1:size(md_ave, 2)]
%    for r = [1:size(fitRange, 1)]
%        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(slope(r, i)));
%    endfor
%endfor

%set(p1, "linewidth", 3);
%set(p2, "linewidth", 3);
%set(p3, "linewidth", 3);
%set(p4, "linewidth", 3);

%title(strcat("Electrical conductivity of 1m NaCl solution - md", "-ave"));
title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",\
             "data interval = ", num2str(skip), ", integration interval = ", num2str(deltaStep)));
%legend("{/Symbol D}t = 1 fs");
xlabel("partial sum upper limit (ps)");
ylabel("Electrical conductivity (S/m)");
axis([0,100,-5, 10]);

print(strcat('ecNoAverageCesaro-', 'ave-slope-skip-', num2str(skip), '-dt-', num2str(deltaStep), '.eps'), '-deps', '-color');
