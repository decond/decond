#!/home/kmtu/bin/octave -qf

clear all;
global numIonTypes;

if (nargin() < 3)
    error("Usage: $plotNoAverageCesaro.m <dataFilename> <skip> <dt>")
endif

dataFilename = argv(){1}
skip = str2num(argv(){2}) #skipped interval in sdCorr data
deltaStep = str2num(argv(){3})

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
%    index = (idx1 - 1) * numIonTypes + idx2;
    # accept only the "upper-half" index pair, because cross-correlation should 
    # be the same for (i,j) and (j,i)
    if (idx1 > idx2)
      error("Error - zipIndexPair: idx1 > idx2");
    else
      index = (idx1 - 1) * numIonTypes + idx2 - (idx1-1);
    endif
endfunction

# md(sdCorr_timeLag, sdCorr_rBin, sdCorrIonTypePairIndex, fileIndex)
load(strcat(dataFilename));
md = D_noAveCesaro(:, :, :);
timeLags = data.timeLags;

numIonTypePairs = size(md, 3) #actually include total part (+1)
num_rBins = size(md, 2)

fitRange = [20, 40; 40, 60; 60, 80; 80, 100]; #ps
%fitRange = [2, 4; 4, 6; 6, 8; 8, 10]; #ps
fitRange *= floor(1000 / skip / deltaStep); #fs (frame)

# calculate slope for each segment of md
for i = [1:numIonTypePairs]
  for j = [1:num_rBins]
    for r = [1:size(fitRange, 1)]
        slope(r,j,i) = polyfit(timeLags(fitRange(r, 1):fitRange(r, 2)), md(fitRange(r, 1):fitRange(r, 2), j, i), 1)(1);
    endfor
  endfor
endfor

save(strcat(dataFilename, '.fit'), "charge", "numIonTypes", "cell", "timestep", "timeLags", "rBins", "rho", "md", "slope");
