#!/home/kmtu/bin/octave -qf

% ***** NOTE *****
% This version2 calculates the md statistics by
% loading the md data files one by one.
% The efficiency is sacrificed to save memory.
% So if the numMD is not large, say, only 10,
% please use the version1, instead.

clear all;
global numIonTypes;

if (nargin() < 4)
    error("Usage: $plotNoAverageCesaro.m <dataFilename> <numMD> <skip> <dt>")
endif

dataFilename = argv(){1}
numMD = str2num(argv(){2})
if (numMD < 2)
    error(cstrcat("Error: numMD < 2. This version2 only accepts multiple MD files.\n",\
                  "       For single MD file, use version1, instead."))
endif

skip = str2num(argv(){3}) #skipped interval in sdCorr data
deltaStep = str2num(argv(){4})

set(0, "defaultlinelinewidth", 4);

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
    index = (idx1 - 1) * numIonTypes + idx2;
endfunction

puts("Loading data files to determine rBins\n");
for n = [1:numMD]
    puts(cstrcat("Loading MD data #", num2str(n), "...\n"));
    dataPath{n} = strcat("./md", num2str(n-1), "/", dataFilename);
    if (n == numMD)
      load(dataPath{n}, "numIonTypes", "timeLags", "rBins");
    else
      load(dataPath{n}, "rBins");
    endif
    num_rBins_tmp(n) = length(rBins);
endfor
num_rBins = min(num_rBins_tmp)
clear("num_rBins_tmp");
rBins = rBins(1:num_rBins);
numIonTypePairs = 1+numIonTypes**2; #actually include total part (+1)

# md(sdCorr_timeLag, sdCorr_rBin, sdCorrIonTypePairIndex, fileIndex)
# calculate md_sum to get md_ave
puts("Loading data files to determine md_sum\n");
md_sum = zeros(length(timeLags), num_rBins, numIonTypePairs);
for n = [1:numMD]
    puts(cstrcat("md_sum: n=", num2str(n), "\n"));
    tmpData = load(dataPath{n}, "ecSDTotalNoAverageCesaro", "ecSDCorrNoAverageCesaro");
    md_sum(:, :, 1) = md_sum(:, :, 1) + tmpData.ecSDTotalNoAverageCesaro(:, 1:num_rBins);
    md_sum(:, :, 2:numIonTypePairs) = md_sum(:, :, 2:numIonTypePairs) + tmpData.ecSDCorrNoAverageCesaro(:, 1:num_rBins, :);
endfor
clear("tmpData");

md_ave = md_sum ./ numMD;
clear("md_sum");

# calculate md_std
md_std = zeros(length(timeLags), num_rBins, numIonTypePairs);
for n = [1:numMD]
    puts(cstrcat("md_std: n=", num2str(n), "\n"));
    tmpData = load(dataPath{n}, "ecSDTotalNoAverageCesaro", "ecSDCorrNoAverageCesaro");
    md_std(:, :, 1) = md_std(:, :, 1) .+ (tmpData.ecSDTotalNoAverageCesaro(:, 1:num_rBins) .- md_ave(:,:,1)).^2;
    md_std(:, :, 2:numIonTypePairs) = md_std(:, :, 2:numIonTypePairs) + (tmpData.ecSDCorrNoAverageCesaro(:, 1:num_rBins, :) .- md_ave(:,:,2:numIonTypePairs)).^2;
endfor
md_std = sqrt(md_std ./ (numMD - 1));
md_err = md_std ./ sqrt(numMD); # standard error


fitRange = [20, 40; 40, 60; 60, 80; 80, 100]; #ps
%fitRange = [2, 4; 4, 6; 6, 8; 8, 10]; #ps
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
slope_b = (S .* Sxy - Sx .* Sy) ./ Delta;
slopeSD = sqrt(S ./ Delta);

save(strcat(dataFilename, '-ave', num2str(numMD), '.fit'), "numIonTypes", "timeLags", "rBins", "md_ave", "md_std", "md_err", "slope", "slopeSD");

%save(strcat('ecNoAverageCesaro-skip-', num2str(skip), '-dt-', num2str(deltaStep), '.fit'), "numIonTypes", "timeLags", "md_ave", "md_std", "md_err", "slope", "slopeSD");

%#numPlots = 1 + numIonTypes + numIonTypes*numIonTypes;
%
%# standard error for selected values
%errFrameInterval = floor(5000 / skip / deltaStep);
%errFrames = [1:errFrameInterval:length(timeLags)]';
%
%
%f1 = figure(1);
%clf;
%hold on;
%errOffset = floor(200 / skip / deltaStep);
%
%for i = [1:numIonTypePairs]
%    drawFrames = errFrames .+ (i-1)*errOffset;
%    drawFrames = drawFrames(drawFrames <= length(timeLags));
%    if (i == 6)
%#        p(i).plot = plot(timeLags, md_ave(:,i), "-", "color", [0, 0.6, 0]);
%        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
%        set(p(i).errorbar(1), "color", [0, 0.6, 0]);
%    elseif (i == 7)
%#        p(i).plot = plot(timeLags, md_ave(:,i), "-", "color", [0, 0, 0]);
%        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i));
%        set(p(i).errorbar(1), "color", [0, 0, 0]);
%    else
%        plotFormat = strcat("-", num2str(i));
%        errorbarFormat = strcat("~", num2str(i));
%        p(i).errorbar = errorbar(timeLags(drawFrames), md_ave(drawFrames,i), md_err(drawFrames,i), errorbarFormat);
%#        p(i).plot = plot(timeLags, md_ave(:,i), plotFormat);
%    endif
%endfor
%
%legend("Total", "Na-Na", "Na-Cl", "Cl-Na", "Cl-Cl", "location", "northwest");
%
%left = 25;
%h_space = 20;
%top = 972;
%v_space = 44;
%for i = [1:size(md_ave, 2)]
%    for r = [1:size(fitRange, 1)]
%        text(left + (r-1)*h_space, top - (i-1)*v_space, num2str(slope(r, i)));
%    endfor
%endfor
%
%title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",\
%             "data interval = ", num2str(skip), ", integration interval = ", num2str(deltaStep)));
%xlabel("partial sum upper limit (ps)");
%ylabel("Non-averaged partial sum (ps*S/m)");
%axis([0,100,-600, 1000]);
%
%print(strcat('ecNoAverageCesaro-', 'ave-skip-', num2str(skip), '-dt-', num2str(deltaStep), '.eps'), '-deps', '-color');
%hold off
%
%##############################
%
%# standard deviation for selected values
%for i = [1:size(fitRange, 1)]
%    stdFrames(i) = (fitRange(i, 1) + fitRange(i, 2)) / 2;
%endfor
%
%f2 = figure(2);
%%clf;
%hold on;
%stdOffset = floor(100 / skip / deltaStep);
%
%for i = [1:numIonTypePairs]
%    drawFrames = stdFrames .+ (i-1)*stdOffset;
%    if (i == 6)
%        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i));
%        set(p(i).errorbar(1), "color", [0, 0.6, 0]);
%    elseif (i == 7)
%        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i));
%        set(p(i).errorbar(1), "color", [0, 0, 0]);
%    else
%        plotFormat = strcat("-", num2str(i));
%        errorbarFormat = strcat("~", num2str(i));
%        p(i).errorbar = errorbar(timeLags(drawFrames), slope(:,i), slopeSD(:,i), errorbarFormat);
%    endif
%endfor
%
%legend("Total", "Auto Na+", "Auto Cl-", "Cross Na-Na", "Cross Na-Cl", "Cross Cl-Na", "Cross Cl-Cl", "location", "northwest");
%
%title(cstrcat("Non-averaged Cesaro sum for electrical conductivity of 1m NaCl solution\n",\
%             "data interval = ", num2str(skip), ", integration interval = ", num2str(deltaStep)));
%xlabel("partial sum upper limit (ps)");
%ylabel("Electrical conductivity (S/m)");
%axis([0,100,-5, 10]);
%
%print(strcat('ecNoAverageCesaro-', 'ave-slope-skip-', num2str(skip), '-dt-', num2str(deltaStep), '.eps'), '-deps', '-color');
