#!/home/kmtu/bin/octave -qf

clear all
format long

if (nargin() < 3)
    error("Usage: $average-sdCorr.m <dataFilename> <numMD> <skip> [num_rBins]");
endif

dataFilename = argv(){1}
numMD = str2num(argv(){2})
skip = str2num(argv(){3}) #skipped interval in sdCorr data

for n = [1:numMD]
    dataPath{n} = strcat("./md", num2str(n-1), "/", dataFilename);
endfor

if (nargin() > 3)
    num_rBins = str2num(argv(){4})
else
    puts("Loading data files to get timeLags, rBins and numIonTypePairs...\n");
    for n = [1:numMD]
        puts(cstrcat("Loading MD data #", num2str(n), "...\n"));
        if (n == numMD)
            load(dataPath{n}, "charge", "timeLags", "rBins");
        else
            load(dataPath{n}, "rBins");
        endif
        num_rBins_tmp(n) = length(rBins);
    endfor
    num_rBins = min(num_rBins_tmp)
    clear("num_rBins_tmp");
endif

# sdCorr(timeLags, rBins, ionTypePairIndex)
# calculate sdCorr_sum to get sdCorr_ave
puts("Loading data files to determine sdCorr_sum\n");
for n = [1:numMD]
    puts(cstrcat("sdCorr_sum: n=", num2str(n), "\n"));
    if (n == 1)
        load(dataPath{n}, "timestep", "charge", "timeLags", "rBins", "sdCorr");
        rBins = rBins(1:num_rBins);
        numIonTypes = length(charge);
%        numIonTypePairs = numIonTypes**2 + 1; #actually include the total part (+1)
        sdCorr_sum = sdCorr(:, 1:num_rBins, :);
%        sdCorr_sum(:, :, 2:numIonTypePairs) = sdCorr(:, 1:num_rBins, :);
%        sdCorr_sum(:, :, 1) = sum(sdCorr(:, 1:num_rBins, :), 3);
    else
        load(dataPath{n}, "sdCorr");
        sdCorr_sum = sdCorr_sum + sdCorr(:, 1:num_rBins, :);
%        sdCorr_sum(:, :, 2:numIonTypePairs) = sdCorr_sum(:, :, 2:numIonTypePairs) + sdCorr(:, 1:num_rBins, :);
%        sdCorr_sum(:, :, 1) = sdCorr_sum(:, :, 1) + sum(sdCorr(:, 1:num_rBins, :), 3);
    endif
endfor
clear("sdCorr");

sdCorr_ave = sdCorr_sum ./ numMD;
clear("sdCorr_sum");

save(strcat(dataFilename, '.ave', num2str(numMD)), "timestep", "charge", "timeLags", "rBins", "sdCorr_ave");

