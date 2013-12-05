#!/home/kmtu/bin/octave -qf

clear all
format long
global numIonTypes

if (nargin() < 2)
    error(cstrcat("Usage: $average-corr.m <data.corr> <numMD>")); 
endif

dataFilename = argv(){1}
numMD = str2num(argv(){2})
%skip = str2num(argv(){3}) #skipped interval in sdCorr data
%extnamePos = rindex(filename, "."); #locate the position of the extension name
%baseFilename = dataFilename(1:extnamePos-1)

load(strcat("./md0/", dataFilename));
clear("autoCorr");
clear("crossCorr");
numIonTypes = length(charge)

function index = zipIndexPair(idx1, idx2)
    global numIonTypes;
    index = (idx1 - 1) * numIonTypes + idx2;
endfunction

function index = zipIndexPair2(idx1, idx2)
    global numIonTypes;
%    index = (idx1 - 1) * numIonTypes + idx2;
    # accept only the "upper-half" index pair, because cross-correlation should 
    # be the same for (i,j) and (j,i)
    if (idx1 > idx2)
      error("Error in zipIndexPair2: idx1 > idx2");
    else
      index = (idx1 - 1) * numIonTypes + idx2 - (idx1-1);
    endif
endfunction

numIonTypePairs = (numIonTypes*(numIonTypes+1)) / 2;

# loading data 
for n = [1:numMD]
    data_tmp = load(strcat("./md", num2str(n-1), "/", dataFilename));
    autoCorrN(:, :, n) = data_tmp.autoCorr;
    crossCorrN(:, :, n) = data_tmp.crossCorr;
endfor
clear("data_tmp");

autoCorr.ave = mean(autoCorrN, 3);
crossCorr.ave = mean(crossCorrN, 3);
autoCorr.std= std(autoCorrN, 0, 3);
crossCorr.std = std(crossCorrN, 0, 3);
autoCorr.err = autoCorr.std / sqrt(numMD);
crossCorr.err = crossCorr.std / sqrt(numMD);

save(strcat(dataFilename, '-ave', num2str(numMD)), "timestep", "charge", "numAtom", "timeLags",\
     "autoCorr", "crossCorr");

