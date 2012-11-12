#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

global totalNumAtoms;

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
        data{i} = readGmx2Matlab_d(vFilename{i});
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
endfor

## data.trajectory(atoms, dimension, frames) 
## squeeze(data{1}.trajectory(:,1,:))' = (frames, atoms) in x-dimension
vCorrTotal = 0;
for dim = [1:3]
    puts(cstrcat("dim = ", num2str(dim), "\n"));
    for n = [1:num_dataFile]
        if (data{n}.num_atoms == 1)
            ## don't transpose because after squeeze it becomes (frames, 1) directly
            vData{n} = squeeze(data{n}.trajectory(:,dim,:)); #(nm / ps)
        else
            vData{n} = squeeze(data{n}.trajectory(:,dim,:))'; #(nm / ps)
        endif
    endfor
    
    puts("calculating vCorrTotal\n");
    vCorrTotal = vCorrTotal + xcorr([vData{:}], maxLag, "unbiased");
endfor
                    
#ex. 3 data files, data{1,2,3}.num_atoms = {2,3,2}
#index = {[1,2],[3,4,5],[6,7]}
#vCorrTotal column: 11,12,...,17,21,22,...,27,...,71,72,...,77
index{1} = [1:data{1}.num_atoms];
for i = [2:num_dataFile]
    index{i} = [index{i-1}(end) + 1: index{i-1}(end) + data{i}.num_atoms];
endfor

function serialIndex = indexPair2SerialIndex(idx1, idx2)
    global totalNumAtoms;
    serialIndex = (idx1 - 1) * totalNumAtoms + idx2;
endfunction


#average 3 dimensions
vCorrTotal = vCorrTotal(maxLag + 1:end, :) / 3;

vAutocorr = cell(1,num_dataFile); #creating cell array
vAutocorr(:) = 0;
vCorr = cell(num_dataFile, num_dataFile);
vCorr(:) = 0;

## loop over all possible index pairs to extract autocorr and corr
for i = [1:num_dataFile]
    for j = [1:num_dataFile]
        for ii = [1:data{i}.num_atoms]
            for jj = [1:data{j}.num_atoms]
                if (i == j && ii == jj)
                    vAutocorr{i} = vAutocorr{i} + vCorrTotal(:,indexPair2SerialIndex(index{i}(ii), index{i}(ii)));
                else
                    vCorr{i,j} = vCorr{i,j} + vCorrTotal(:,indexPair2SerialIndex(index{i}(ii), index{j}(jj)));
                endif
            endfor
        endfor
    endfor
endfor

# we need number of atoms for each ion type to calculate diffusion coefficients
for i = [1:num_dataFile]
    numAtoms(i) = data{i}.num_atoms;
endfor

# output time vector for convenience of plotting
timeLags = [0:maxLag]' * timestep;

save(strcat(outFilename, ".vCorr"), "timestep", "charge", "numAtoms", "timeLags", "vAutocorr", "vCorr");

vCorrTotal = sum(vCorrTotal, 2);
save(strcat(outFilename, ".vCorrTotal"), "vCorrTotal");

