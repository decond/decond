#!/home/kmtu/bin/octave -qf

clear all
format long
global numIonTypes

if (nargin() < 2)
    error(cstrcat("Usage: $average-vCorr.m <data.vCorr> <numMD>")); 
endif

dataFilename = argv(){1}
numMD = str2num(argv(){2})
%skip = str2num(argv(){3}) #skipped interval in sdCorr data
%extnamePos = rindex(filename, "."); #locate the position of the extension name
%baseFilename = dataFilename(1:extnamePos-1)


for n = [1:numMD]
    dataPath{n} = strcat("./md", num2str(n-1), "/", dataFilename);
    md(n) = load(dataPath{n});
endfor

timestep = md(1).timestep;
charge = md(1).charge;
numAtoms = md(1).numAtoms;
timeLags = md(1).timeLags;
numIonTypes = length(charge)

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

vAutocorr_ave = zeros(length(timeLags), numIonTypes);
vCorr_ave = zeros(length(timeLags), numIonTypes*(numIonTypes+1)/2);
%vAutocorr = cell(size(md(1).vAutocorr));
%vCorr = cell(size(md(2).vCorr));
%vAutocorr(:) = 0;
%vCorr(:) = 0;

% sum
for n = [1:numMD]
    for i = [1:numIonTypes]
        vAutocorr_ave(:,i) += md(n).vAutocorr{i};
        for j = [i:numIonTypes]
            idx = zipIndexPair(i,j);
            vCorr_ave(:,idx) .+= (md(n).vCorr{i,j} .+ md(n).vCorr{j,i}) ./ 2;
        endfor
    endfor
endfor

% average    
vAutocorr_ave /= numMD;
vCorr_ave /= numMD;
%for i = [1:length(vAutocorr)]
%    vAutocorr{i} = vAutocorr{i} / numMD;
%    for j = [1:length(vAutocorr)]
%        vCorr{i,j} = vCorr{i,j} / numMD;
%    endfor
%endfor

vAutocorr_std = zeros(size(vAutocorr_ave));
vCorr_std = zeros(size(vCorr_ave));
%std
for n = [1:numMD]
    for i = [1:numIonTypes]
        vAutocorr_std(:,i) += (md(n).vAutocorr{i} .- vAutocorr_ave(:, i)).**2;
        for j = [i:numIonTypes]
            idx = zipIndexPair(i,j);
            vCorr_std(:,idx) += ((md(n).vCorr{i,j} + md(n).vCorr{j,i})/2 .- vCorr_ave(:,idx)).**2;
        endfor
    endfor
endfor
vAutocorr_std /= (numMD - 1);
vCorr_std /= (numMD - 1);
vAutocorr_std = sqrt(vAutocorr_std);
vCorr_std = sqrt(vCorr_std);

vAutocorr_err = vAutocorr_std ./ sqrt(numMD);
vCorr_err = vCorr_std ./ sqrt(numMD);

save(strcat(dataFilename, '-ave', num2str(numMD)), "timestep", "charge", "numAtoms", "timeLags",\
     "vAutocorr_ave", "vCorr_ave", "vAutocorr_err", "vCorr_err" );

