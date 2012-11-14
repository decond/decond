#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

global totalNumAtoms dataIndex a2g;

num_parArg = 2;
num_dataArg = nargin() - num_parArg;
if (num_dataArg < 2 || mod(num_dataArg, 2) == 1)
    error("Usage: $calculateCorr.m <outFilename> <maxLag (-1=max)> <velData1.binary> <charge1> [<velData2.binary> <charge2>...]")
else
    outFilename = argv(){1};

    maxLag = str2num(argv(){2}); # in the unit of frame number

    num_dataFile = num_dataArg / 2
    for i = [1: num_dataFile]
        vFilename{i} = argv(){num_parArg + 2*i - 1};
        charge(i) = str2num(argv(){num_parArg + 2*i});
        data{i} = readGmx2Matlab(vFilename{i});
    endfor
endif

## check the num_frames are the same for all data
for n = [1:num_dataFile-1]
    if (data{n}.num_frames != data{n+1}.num_frames)
        error(cstrcat("Numbers of frames are different between ", vFilename{n}, " and ", vFilename{n+1}))
    endif
    if (data{n}.time_step != data{n+1}.time_step)
        error(cstrcat("Timesteps are different between ", vFilename{n}, " and ", vFilename{n+1}))
    endif
endfor

timestep = data{1}.time_step
num_frames = data{1}.num_frames #for showing purpose
if (maxLag < 0)
    maxLag = num_frames - 1;
endif
maxLag #showing

totalNumAtoms = 0;
for n = [1:num_dataFile]
    totalNumAtoms = totalNumAtoms + data{n}.num_atoms; 
    # we need number of atoms for each ion type to calculate diffusion coefficients
    numAtoms(n) = data{n}.num_atoms;
endfor

#ex. 3 data files, data{1,2,3}.num_atoms = {2,3,2}
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

    if (!exist("data"))
        for i = [1: num_dataFile]
            data{i} = readGmx2Matlab(vFilename{i});
        endfor
    endif

    for n = [1:num_dataFile]
        ## data.trajectory(atoms, dimension, frames) 
        ## squeeze(data{1}.trajectory(:,1,:))' = (frames, atoms) in x-dimension
        if (numAtoms(n) == 1)
            ## don't transpose because after squeeze it becomes (frames, 1) directly
            vData{n} = squeeze(data{n}.trajectory(:,dim,:)); #(nm / ps)
        else
            vData{n} = squeeze(data{n}.trajectory(:,dim,:))'; #(nm / ps)
        endif
    endfor

    clear data;
    
    puts("calculating correlation\n");
    whos
    [vAutocorr_tmp, vCorr_tmp] =  xcorr_tu([vData{:}], maxLag, "unbiased");
    for i = [1:num_dataFile]
            vAutocorr{i} = vAutocorr{i} + vAutocorr_tmp{i};
            vCorrTotal = vCorrTotal + vAutocorr_tmp{i};
        for j = [1:num_dataFile]
            vCorr{i,j} = vCorr{i,j} + vCorr_tmp{i,j};
            vCorrTotal = vCorrTotal + vCorr_tmp{i,j};
        endfor
    endfor
endfor
                    
#average 3 dimensions
#vCorrTotal = vCorrTotal(maxLag + 1:end, :) / 3;
vCorrTotal = vCorrTotal / 3;
for i = [1:num_dataFile]
    vAutocorr{i} = vAutocorr{i} / 3;
    for j = [1:num_dataFile]
        vCorr{i,j} = vCorr{i,j} / 3;
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

