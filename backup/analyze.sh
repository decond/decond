if [[ ! $# == 1 ]]; then
  echo "Usage: $0 <md index>"
  exit 1;
fi

idx=$1
echo "md index = $idx"

if [[ ! -d "md$idx" ]]; then
  echo "md$idx is not found"
  exit 1;
fi

cd md$idx
boxLength=$(tail -n 1 frame$idx.gro | awk '{print $1}')
systemVolume=$(tail -n 1 frame$idx.gro | echo "$(awk '{print $1}')^3" | bc)
echo "box length = $boxLength"
echo "system volume = $systemVolume"

GMDIR=/home/kmtu/local/gromacs-4.5.5/bin
#NDXDIR=/home/kmtu/elec_cond/NaCl/s36-w1000/NVE
ECDIR=/home/kmtu/elec_cond/econd

mdData=md.trr
naBase=md-na
clBase=md-cl
naData=$naBase.trr
clData=$clBase.trr
naXBin=${naBase}_d_xdata.binary
clXBin=${clBase}_d_xdata.binary
naVBin=${naBase}_d_vdata.binary
clVBin=${clBase}_d_vdata.binary
dataMaxlag=1000000
outName=lag$dataMaxlag
analyzeMaxlag=100000

#$GMDIR/trjconv_d -f $mdData -o $naData -n ../na.ndx &&\
#$GMDIR/trjconv_d -f $mdData -o $clData -n ../cl.ndx &&\
#$ECDIR/convertData.m $naData $naBase double &&\
#$ECDIR/convertData.m $clData $clBase double &&\
#$ECDIR/calCorr.m $outName $dataMaxlag double $naVBin 1 $clVBin -1 &&\
$ECDIR/spatialDecompose.m $outName $dataMaxlag 0.1 $boxLength double $naXBin $naVBin 1 $clXBin $clVBin -1
#$ECDIR/calDC.m $outName.vCorr $analyzeMaxlag &&\
#$ECDIR/calEC.m $outName.vCorr $analyzeMaxlag $systemVolume &&\
#mv $outName.dc ${outName}-int_$analyzeMaxlag.dc &&\
#mv $outName.ec ${outName}-int_$analyzeMaxlag.ec 
#g_rdf_d -f $mdData -n ../na-cl.ndx -o rdf-na_cl.xvg
#$ECDIR/calNoAverageCesaroEC.m $outName.vCorr $analyzeMaxlag $systemVolume

cd ..
