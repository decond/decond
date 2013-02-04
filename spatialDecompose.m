#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

global totalNumAtoms dataIndex a2g;

num_parArg = 4;
num_argPerData = 3;
num_dataArg = nargin() - num_parArg;
if (num_dataArg < num_argPerData || mod(num_dataArg, num_argPerData) != 0)
    error("Usage: $spatialDecompose.m <outFilename> <maxLag (-1=max)> <rBinWidth(nm)> <precision> <posData1.binary> <velData1.binary> <charge1> [<posData2.binary> <velData2.binary> <charge2>...]")
else
    outFilename = argv(){1};
    maxLag = str2num(argv(){2}); # in the unit of frame number
    rBinWidth = str2num(argv(){3}); 
    precision = argv(){4}; # single or double

    num_dataFile = num_dataArg / num_argPerData
    for i = [1: num_dataFile]
        rFilename{i} = argv(){num_parArg + num_argPerData*(i-1) + 1};
        vFilename{i} = argv(){num_parArg + num_argPerData*(i-1) + 2};
        charge(i) = str2num(argv(){num_parArg + num_argPerData*(i-1) + 3});
        rData{i} = readGmx2Matlab_tu(rFilename{i}, precision);
        vData{i} = readGmx2Matlab_tu(vFilename{i}, precision);
    endfor
endif

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
for i = [2:num_dataFile]
    dataIndex{i} = [dataIndex{i-1}(end) + 1: dataIndex{i-1}(end) + numAtoms(i)];
endfor

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

a2g = @atomIndex2GroupIndex;

vCorrTotal = 0;
vAutocorr = cell(1,num_dataFile); #creating cell array
vAutocorr(:) = 0;
vCorr = cell(num_dataFile, num_dataFile);
vCorr(:) = 0;

for dim = [1:3]
    puts(cstrcat("dim = ", num2str(dim), "\n"));

    #Not sure why the check is here, is it necessary?
    if (!exist("vData"))
        for i = [1: num_dataFile]
            vData{i} = readGmx2Matlab_tu(vFilename{i}, precision);
        endfor
    endif

    for n = [1:num_dataFile]
        ## data.trajectory(atoms, dimension, frames) 
        ## squeeze(vData{1}.trajectory(:,1,:))' = (frames, atoms) in x-dimension
        if (numAtoms(n) == 1)
            ## don't transpose because after squeeze it becomes (frames, 1) directly
            vDataSqz{n} = squeeze(vData{n}.trajectory(:,dim,:)); #(nm / ps)
        else
            vDataSqz{n} = squeeze(vData{n}.trajectory(:,dim,:))'; #(nm / ps)
        endif
    endfor

    clear data;
    
    puts("calculating correlation\n");
    whos
    [vAutocorr_tmp, vCorr_tmp] =  xcorr_tu([vDataSqz{:}], maxLag, "unbiased");
    for i = [1:num_dataFile]
            vAutocorr{i} = vAutocorr{i} + vAutocorr_tmp{i};
#            vCorrTotal = vCorrTotal + vAutocorr{i};
        for j = [1:num_dataFile]
            vCorr{i,j} = vCorr{i,j} + vCorr_tmp{i,j};
#            vCorrTotal = vCorrTotal + vCorr_tmp{i,j};
        endfor
    endfor
endfor

#average 3 dimensions
#vCorrTotal = vCorrTotal(maxLag + 1:end, :) / 3;
#vCorrTotal = vCorrTotal / 3;
for i = [1:num_dataFile]
    vAutocorr{i} = vAutocorr{i} / 3;
    vCorrTotal = vCorrTotal + vAutocorr{i};
    for j = [1:num_dataFile]
        vCorr{i,j} = vCorr{i,j} / 3;
        vCorrTotal = vCorrTotal + vCorr{i,j};
    endfor
endfor

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

save(strcat(outFilename, ".vCorr"), "timestep", "charge", "numAtoms", "timeLags", "vAutocorr", "vCorr");

#vCorrTotal = sum(vCorrTotal, 2);
save(strcat(outFilename, ".vCorrTotal"), "vCorrTotal");

