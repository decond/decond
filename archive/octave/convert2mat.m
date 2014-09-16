#!/home/kmtu/bin/octave -qf

clear all
format long

if (nargin() != 2)
    error("Usage: $convert2mat.m <octData> <outFileName>")
endif

octData = argv(){1}
outFileName = argv(){2}

load(octData)
save('-7', outFileName)

