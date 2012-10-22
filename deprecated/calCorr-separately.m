#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

num_parArg = 2;
num_dataArg = nargin() - num_parArg;
if (num_dataArg < 2 || mod(num_dataArg, 2) == 1)
    error("Usage: $calculateCorr.m <outFilename> <maxLag (-1=max)> <velData1.binary> <charge1> [<velData2.binary> <charge2>...]")
else
    outFilename = argv(){1};

    maxLag = str2num(argv(){2}) # in the unit of frame number

    num_dataFile = num_dataArg / 2
    for i = [1: num_dataFile]
        vFilename{i} = argv(){num_parArg + 2*i - 1};
        charge{i} = str2num(argv(){num_parArg + 2*i});
        data{i} = readGmx2Matlab(vFilename{i});
    endfor
endif

#check the num_frames are the same for all data
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
    maxLag = num_frames - 1
endif

## data.trajectory(atoms, dimension, frames) 
## squeeze(data{1}.trajectory(:,1,:))' = (frames, atoms) in x-dimension
jAutocorr = cell(1,num_dataFile); #creating cell array
jAutocorr(:) = 0;
jCorr = cell(num_dataFile, num_dataFile);
jCorr(:) = 0;
for dim = [1:3]
    puts(cstrcat("dim = ", num2str(dim), "\n"));
    for n = [1:num_dataFile]
        if (data{n}.num_atoms == 1)
            #don't transpose because after squeeze it becomes (frames, 1) directly
            jData{n} = charge{n}*squeeze(data{n}.trajectory(:,dim,:)); #(nm / ps)
        else
            jData{n} = charge{n}*squeeze(data{n}.trajectory(:,dim,:))'; #(nm / ps)
        endif
    endfor
    
#For test purpose
#    #calculate total diffusion directly
#    puts("calculating jCorrTotal\n");
#    if (dim == 1)
#        [jCorrTotal, t] = xcorr([jData{:}], maxLag, "unbiased");
#    else
#        jCorrTotal = jCorrTotal + xcorr([jData{:}], maxLag, "unbiased");
#    endif
#################


    #calculate each term separately
    for i = [1:num_dataFile]
        for j = [1:num_dataFile]
            for ii = [1:data{i}.num_atoms]
                if (i == j)
                    #self-diffusion (jAutocorrelation) for data{i}
                    puts(cstrcat("calculating jAutocorr", num2str(i), " of dimension ", num2str(dim), " for atom ", num2str(ii), "\n"));
                    jAutocorr{i} = jAutocorr{i} + xcorr(jData{i}(:, ii), maxLag, "unbiased");
                endif
                for jj = [1:data{j}.num_atoms]
                    if (i == j)
                        if (ii != jj)
                            #mutual-diffusion for data{i}-data{i}
                            puts(cstrcat("calculating jCorr", num2str(i), num2str(i), " of dimension ", num2str(dim), " for atom-pair(", num2str(ii), ", ", num2str(jj), ")\n"));
                            jCorr{i,j} = jCorr{i,j} + xcorr(jData{i}(:,ii), jData{j}(:,jj), maxLag, "unbiased");
                        endif
                    else
                            #mutual-diffusion for data{i}-data{j}
                            puts(cstrcat("calculating jCorr", num2str(i), num2str(i), " of dimension ", num2str(dim), " for atom-pair(", num2str(ii), ", ", num2str(jj), ")\n"));
                            jCorr{i,j} = jCorr{i,j} + xcorr(jData{i}(:,ii), jData{j}(:,jj), maxLag, "unbiased");
                    endif
                endfor
            endfor
        endfor
    endfor
endfor
                    

#average 3 dimensions
#jCorrTotal = sum(jCorrTotal(maxLag + 1:end, :), 2) / 3;
#save(strcat(outFilename, ".jCorrTotal"), "jCorrTotal");

#jCorrSum = 0;
for i = [1:num_dataFile]
    jAutocorr{i} = jAutocorr{i}(maxLag + 1: end) / 3;
#    jCorrSum = jCorrSum + jAutocorr{i};
    for j = [1:num_dataFile]
        if (length(jCorr{i,j}) > 1 || jCorr{i,j} > 0)
            jCorr{i,j} = jCorr{i,j}(maxLag + 1: end) / 3;
#            jCorrSum = jCorrSum + jCorr{i, j}; #For test
        endif
    endfor
endfor
#save(strcat(outFilename, ".jAutocorr"), "jAutocorr");
#save(strcat(outFilename, ".jCorr"), "jCorr");
save(strcat(outFilename, ".jCorr"), "timestep", "jAutocorr", "jCorr");

#For test
#save(strcat(outFilename, ".jCorrSum"), "jCorrSum");

