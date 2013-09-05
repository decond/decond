if [[ ! $# == 2 ]]; then
  echo "Usage: $0 <md index> <trajInterval>"
  exit 1;
fi

idx=$1
echo "md index = $idx"

if [[ ! -d "md$idx" ]]; then
  echo "md$idx is not found"
  exit 1;
fi

trajInterval=$2
echo "trajInterval = $trajInterval"

cd md$idx
boxLength=$(tail -n 1 frame$idx.gro | awk '{print $1}')
systemVolume=$(tail -n 1 frame$idx.gro | echo "$(awk '{print $1}')^3" | bc)
echo "box length = $boxLength"
echo "system volume = $systemVolume"

#GMDIR=/home/kmtu/local/gromacs-4.5.5/bin
GMDIR=/home/kmtu/local/gromacs-4.6.3/bin
#NDXDIR=/home/kmtu/elec_cond/NaCl/s36-w1000/NVE
ECDIR=/home/kmtu/elec_cond/econd
TRJCONV=$GMDIR/trjconv_mpi_d

mdData=md-na-cl.trr
naBase=md-na-skip$trajInterval
clBase=md-cl-skip$trajInterval
naData=$naBase.trr
clData=$clBase.trr
naXBin=${naBase}_d_xdata.binary
clXBin=${clBase}_d_xdata.binary
naVBin=${naBase}_d_vdata.binary
clVBin=${clBase}_d_vdata.binary

maxlag=1000000
dataMaxlag=$((maxlag/trajInterval))
#outName=lag$maxlag-skip$trajInterval
outName=E6E5-skip$trajInterval

analyzeMaxlag=$((100000/trajInterval))

#$TRJCONV -f md-na-cl.g96 -o $mdData
#$TRJCONV -f $mdData -o $naData -n ../na-pure.ndx -skip $trajInterval &&\
#$TRJCONV -f $mdData -o $clData -n ../cl-pure.ndx -skip $trajInterval 

#$ECDIR/convertData.m $naData $naBase double &&\
#$ECDIR/convertData.m $clData $clBase double 
#$ECDIR/calCorr.m $outName $dataMaxlag double $naVBin 1 $clVBin -1

#$ECDIR/spatialDecompose.m $outName $dataMaxlag 0.1 $boxLength double $naXBin $naVBin 1 $clXBin $clVBin -1
#$ECDIR/calDC.m $outName.vCorr $analyzeMaxlag &&\
#$ECDIR/calEC.m $outName.vCorr $analyzeMaxlag $systemVolume &&\
#mv $outName.dc ${outName}-int_$analyzeMaxlag.dc &&\
#mv $outName.ec ${outName}-int_$analyzeMaxlag.ec 
#g_rdf_d -f $mdData -n ../na-cl.ndx -o rdf-na_cl.xvg

#for d in 15 20 25 30 35 40 
#do
    d=1
    $ECDIR/calNoAverageCesaroEC.m $outName.vCorr $analyzeMaxlag $systemVolume $d
#done

cd ..
