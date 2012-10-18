#!/home/kmtu/bin/octave -qf

clear all
format long
pkg load signal;

num_parArg = 2;
num_dataArg = nargin() - num_parArg;
if (num_dataArg < 2 || mod(num_dataArg, 2) == 1)
    error("Usage: $calculateCorr.m <outFileName> <maxLag (-1=max)> <velData1.binary> <charge1> [<velData2.binary> <charge2>...]")
else
    outFileName = argv(){1};

    maxLag = str2num(argv(){2}) # in the unit of frame number

    num_dataFile = num_dataArg / 2
    for i = [1: num_dataFile]
        vFileName{i} = argv(){num_parArg + 2*i - 1};
        charge{i} = str2num(argv(){num_parArg + 2*i});
        data{i} = readGmx2Matlab(vFileName{i});
    endfor
endif

#check the num_frames are the same for all data
for n = [1:num_dataFile-1]
    if (data{n}.num_frames != data{n+1}.num_frames)
        error(cstrcat("Numbers of frames are different between ", vFileName{n}, " and ", vFileName{n+1}))
    endif
endfor

if (maxLag < 0)
    maxLag = data{1}.num_frames - 1
endif

## data.trajectory(atoms, dimension, frames) 
## squeeze(data{1}.trajectory(:,1,:))' = (frames, atoms) in x-dimension
autocorrData = cell(1,num_dataFile); #creating cell array
autocorrData(:) = 0;
corrData = cell(num_dataFile, num_dataFile);
corrData(:) = 0;
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
    #calculate total diffusion directly
    puts("calculating corrData\n");
    if (dim == 1)
        [corrDataTotal, t] = xcorr([jData{:}], maxLag, "unbiased");
    else
        corrDataTotal = corrDataTotal + xcorr([jData{:}], maxLag, "unbiased");
    endif
#################


    #calculate each term separately
    for i = [1:num_dataFile]
        for j = [1:num_dataFile]
            for ii = [1:data{i}.num_atoms]
                if (i == j)
                    #self-diffusion (autocorrelation) for data{i}
                    puts(cstrcat("calculating autocorrData", num2str(i), " of dimension ", num2str(dim), " for atom ", num2str(ii), "\n"));
                    autocorrData{i} = autocorrData{i} + xcorr(jData{i}(:, ii), maxLag, "unbiased");
                endif
                for jj = [1:data{j}.num_atoms]
                    if (i == j)
                        if (ii != jj)
                            #mutual-diffusion for data{i}-data{i}
                            puts(cstrcat("calculating corrData", num2str(i), num2str(i), " of dimension ", num2str(dim), " for atom-pair(", num2str(ii), ", ", num2str(jj), ")\n"));
                            corrData{i,j} = corrData{i,j} + xcorr(jData{i}(:,ii), jData{j}(:,jj), maxLag, "unbiased");
                        endif
                    else
                            #mutual-diffusion for data{i}-data{j}
                            puts(cstrcat("calculating corrData", num2str(i), num2str(i), " of dimension ", num2str(dim), " for atom-pair(", num2str(ii), ", ", num2str(jj), ")\n"));
                            corrData{i,j} = corrData{i,j} + xcorr(jData{i}(:,ii), jData{j}(:,jj), maxLag, "unbiased");
                    endif
                endfor
            endfor
        endfor
    endfor
endfor
                    

#average 3 dimensions
corrDataTotal = corrDataTotal(maxLag + 1:end) / 3;
save(strcat(outFileName, ".corrDataTotal"), "corrDataTotal");

sumCorrData = 0
for i = [1:num_dataFile]
    autocorrData{i} = autocorrData{i}(maxLag + 1: end) / 3;
    sumCorrData = sumCorrData + autocorrData{i}
    for j = [1:num_dataFile]
        corrData{i,j} = corrData{i,j}(maxLag + 1: end) / 3;
        sumCorrData = sumCorrData + corrData{i} #For test
    endfor
endfor
save(strcat(outFileName, ".autocorrData"), "autocorrData");
save(strcat(outFileName, ".corrData"), "corrData");

#For test
save(strcat(outFileName, ".sumCorrData"), "sumCorrData");

