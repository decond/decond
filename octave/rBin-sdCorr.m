#!/home/kmtu/bin/octave -qf

clear all
format long

if (nargin() < 2)
    error(cstrcat("Usage: $rBin-sdCorr.m <data.sdcorr> <rBinWindow>\n",\
                   "rBinWindow: combine every <rBinWindow> into one bin"));
endif

dataFilename = argv(){1}
rBinWindow = str2num(argv(){2})

puts(cstrcat("Loading sdcorr data ...\n"));
load(dataFilename, "timestep", "charge", "numAtom", "volume_ave", "timeLags", "rBins", "sdCorr_ave", "rho_ave");

numIonTypes = length(charge);
numIonTypePairs = numIonTypes**2;

rBins_source = rBins;
clear("rBins");

num_rBins = floor(length(rBins_source) / rBinWindow);
rBins = zeros(num_rBins, 1);
sdCorr = zeros(size(sdCorr_ave, 1), num_rBins, size(sdCorr_ave, 3));
rho = zeros(num_rBins, size(rho_ave, 2));

for i = [1 : rBinWindow]
  idx_range = [i : rBinWindow : num_rBins * rBinWindow];
  rBins += rBins_source(idx_range);
  sdCorr += sdCorr_ave(:, idx_range, :);
  rho += rho_ave(idx_range, :);
endfor
rBins /= rBinWindow;
sdCorr /= rBinWindow;
rho /= rBinWindow;

save(strcat(dataFilename, '-rBinWindow', num2str(rBinWindow)), "timestep", "charge",\
     "volume_ave", "numAtom", "timeLags", "rBins", "sdCorr", "rho");
