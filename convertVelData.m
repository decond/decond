#!/home/kmtu/bin/octave -qf

clear all
format long

if (nargin() < 3)
    error("Usage: $convertVelData.m <data.trr> <outFileName> <precision>")
endif

dataFileName = argv(){1}
outFileName = argv(){2}
fileType = argv(){3}

trr2matlab_tu(dataFileName, 'v', outFileName, fileType)
#trr2matlab(dataFileName, 'v', outFileName)

