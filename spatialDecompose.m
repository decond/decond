#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

global totalNumAtoms dataIndex a2g rBinWidth boxLength;

num_parArg = 5;
num_argPerData = 3;
num_dataArg = nargin() - num_parArg;
if (num_dataArg < num_argPerData || mod(num_dataArg, num_argPerData) != 0)
    error("Usage: $spatialDecompose.m <outFilename> <maxLag (-1=max)> <rBinWidth(nm)> <boxLength(nm)> <precision> <posData1.binary> <velData1.binary> <charge1> [<posData2.binary> <velData2.binary> <charge2>...]")
else
    outFilename = argv(){1};
    maxLag = str2num(argv(){2}); # in the unit of frame number
    rBinWidth = str2num(argv(){3}); 
    boxLength = str2num(argv(){4}); 
    precision = argv(){5}; # single or double

    num_rBin = ceil(boxLength / rBinWidth)

    num_dataFile = num_dataArg / num_argPerData
    for i = [1: num_dataFile]
        rFilename{i} = argv(){num_parArg + num_argPerData*(i-1) + 1};
        vFilename{i} = argv(){num_parArg + num_argPerData*(i-1) + 2};
        charge(i) = str2num(argv(){num_parArg + num_argPerData*(i-1) + 3});
        rData{i} = readGmx2Matlab_tu(rFilename{i}, precision);
        vData{i} = readGmx2Matlab_tu(vFilename{i}, precision);
    endfor
endif

puts("Tag1\n");
whos

## check the num_frames are the same for all data
for n = [1:num_dataFile-1]
    if (rData{n}.num_frames != rData{n+1}.num_frames)
        error(cstrcat("Numbers of frames are different between ", rFilename{n}, " and ", rFilename{n+1}))
    endif
    if (rData{n}.time_step != rData{n+1}.time_step)
        error(cstrcat("Timesteps are different between ", rFilename{n}, " and ", rFilename{n+1}))
    endif
    if (vData{n}.num_frames != vData{n+1}.num_frames)
        error(cstrcat("Numbers of frames are different between ", vFilename{n}, " and ", vFilename{n+1}))
    endif
    if (vData{n}.time_step != vData{n+1}.time_step)
        error(cstrcat("Timesteps are different between ", vFilename{n}, " and ", vFilename{n+1}))
    endif
endfor

timestep = vData{1}.time_step
num_frames = vData{1}.num_frames #for showing purpose
if (maxLag < 0)
    maxLag = num_frames - 1;
endif
maxLag #showing

totalNumAtoms = 0;
for n = [1:num_dataFile]
    totalNumAtoms = totalNumAtoms + vData{n}.num_atoms; 
    # we need number of atoms for each ion type to calculate diffusion coefficients
    numAtoms(n) = vData{n}.num_atoms;
endfor


#ex. 3 data files, vData{1,2,3}.num_atoms = {2,3,2}
#dataIndex = {[1,2],[3,4,5],[6,7]}
#atomIndex = 1,2,3,4,5,6,7
#groupIndex = [1], [2], [3]
#---- after doing correlation ----
#vCorrTotal column (atomIndex pair): 11,12,...,17,21,22,...,27,...,71,72,...,77
#serialIndex: 1, 2,..., 7, 8, 9,... 49
#groupIndex pair: [1][1]=11,12,21,22=(1,2)x(1,2); [1][2]=(1,2)x(3,4,5); ...

dataIndex{1} = [1:numAtoms(1)];
for n = [2:num_dataFile]
    dataIndex{n} = [dataIndex{n-1}(end) + 1: dataIndex{n-1}(end) + numAtoms(n)];
endfor

#rData{i} should have only one time frame. combine all ion types to make one array
rDataAll = [];
for n = [1:num_dataFile]
    rDataAll = [rDataAll; rData{n}.trajectory(:,:,1)];
endfor
clear rData;
if (size(rDataAll,1) != totalNumAtoms)
    error("Combining rData to rDataAll failed. The total numbers of atoms are inconsistent");
endif

puts("Tag2\n");
whos

function serialIndex = atomPair2SerialIndex(idx1, idx2)
    global totalNumAtoms;
    serialIndex = (idx1 - 1) * totalNumAtoms + idx2;
endfunction

function groupIndex = atomIndex2GroupIndex(idx)
    global dataIndex;
    for i = [1:length(dataIndex)]
        if (any(dataIndex{i} == idx))
            groupIndex = i;
            return;
        endif
    endfor
    error(strcat("Unable to convert atomIndex:", num2str(idx), " to groupIndex"));
endfunction

function wrappedR = wrap(r)
    global boxLength;
    wrappedR = r;
    isOutsideHalfBox = (r > (boxLength / 2));
    wrappedR(isOutsideHalfBox) = abs(r(isOutsideHalfBox) - boxLength);
endfunction

function rBinIndex = getBinIndex(pos1, pos2)
    global rBinWidth;
    r = abs(pos1 - pos2);
    r = wrap(r);
    r = sqrt(sum(r.*r, 2));
    rBinIndex = ceil(r ./ rBinWidth);
    rBinIndex(rBinIndex == 0) = 1;
endfunction

a2g = @atomIndex2GroupIndex;

#cAutocorr = cell(1,num_dataFile); #creating cell array
#cAutocorr(:) = zeros(maxLag+1, num_rBin);
cCorr = cell(num_dataFile, num_dataFile);
cCorr(:) = zeros(maxLag+1, num_rBin);
#rhoAutocorr = cell(1, num_dataFile); #rho_I(r)
#rhoAutocorr(:) = zeros(num_rBin);
rhoCorr = cell(num_dataFile, num_dataFile); #rho_IJ(r)
rhoCorr(:) = zeros(1, num_rBin);

#rBinIndex(i,j) = the binning index of |r_i - r_j|, the distance between particle i and j
rBinIndex = zeros(totalNumAtoms);
for k = [0:totalNumAtoms-1]
    rBinIndexK = getBinIndex(rDataAll(1+k:end,:), rDataAll(1:end-k,:));
    for i = [1:length(rBinIndexK)]
        rBinIndex(i, i+k) = rBinIndexK(i);
    endfor
endfor
rBinIndex = rBinIndex + triu(rBinIndex, 1)' #mirror the upper half to lower half
clear rDataAll;

puts("Tag3\n");
whos

## vData{i}.trajectory(atoms, dimension, frames) 
for i = [1:num_dataFile]
    for j = [1:num_dataFile]
        for ii = [1:vData{i}.num_atoms]
#            if (i == j)
#                idx = rBinIndex(dataIndex{i}(ii), dataIndex{i}(ii));
#                for dim = [1:3]
#                    cAutocorr{i}(:, idx) += xcorr(squeeze(vData{i}.trajectory(ii, dim, :)), maxLag, "unbiased")(maxLag+1:end);
#                    rhoAutocorr{i}(idx) += 1;
#                endfor 
#            endif
            for jj = [1:vData{j}.num_atoms]
                #omit autocorrelation
                if ((i != j) || (ii != jj))
                    for dim = [1:3]
                        idx = rBinIndex(dataIndex{i}(ii), dataIndex{j}(jj));
                        cCorr{i,j}(:, idx) += xcorr(squeeze(vData{i}.trajectory(ii, dim, :)), squeeze(vData{j}.trajectory(jj, dim, :)), maxLag, "unbiased")(maxLag+1:end);
                        rhoCorr{i,j}(idx) += 1;
                    endfor
                endif
            endfor 
        endfor
    endfor
endfor

clear vData;
puts("Tag4\n");
whos

#average 3 dimensions
for i = [1:num_dataFile]
#    cAutocorr{i} ./= (3*rhoAutocorr{i});
    for j = [1:num_dataFile]
        cCorr{i,j} ./= (3*rhoCorr{i,j});
    endfor
endfor

puts("Tag5\n");
whos

### loop over all possible dataIndex pairs to extract autocorr and corr
#for i = [1:num_dataFile]
#    for j = [1:num_dataFile]
#        for ii = [1:numAtoms(i)]
#            for jj = [1:numAtoms(j)]
#                if (i == j && ii == jj)
#                    vAutocorr{i} = vAutocorr{i} + vCorrTotal(:,atomPair2SerialIndex(dataIndex{i}(ii), dataIndex{i}(ii)));
#                else
#                    vCorr{i,j} = vCorr{i,j} + vCorrTotal(:,atomPair2SerialIndex(dataIndex{i}(ii), dataIndex{j}(jj)));
#                endif
#            endfor
#        endfor
#    endfor
#endfor


# output time vector for convenience of plotting
timeLags = [0:maxLag]' * timestep;
rBins = ([1:num_rBin] - 0.5)* rBinWidth;

save(strcat(outFilename, ".cCorr"), "timestep", "charge", "numAtoms", "timeLags", "rBins", "cCorr", "rhoCorr");
#save(strcat(outFilename, ".cCorrUnaverage"), "timestep", "charge", "numAtoms", "timeLags", "rBins", "cCorr", "rhoCorr");

