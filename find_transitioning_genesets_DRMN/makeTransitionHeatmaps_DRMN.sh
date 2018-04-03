#!/bin/bash
# Modified for DRMN to get appropriate colormaps for expression data  

if [ $# -lt 5 ]
then
	echo "Too few arguments $#"
	echo "Usage: ./makeTransitionHeatmaps.sh outputdir numberofclusters minExpr maxExpr heatmap_script_path"
	exit
fi
OUTDIR=$1
FIGDIR=$OUTDIR/figs
if [ -d $FIGDIR ]
then
	echo "Please remove the figs directory in $OUTDIR"
	exit
fi
mkdir $FIGDIR
#MAXCID=`expr $2 - 1` # actually, these are getting offset by 1
MAXCID=$2
HALFCID=`echo $MAXCID| awk '{printf("%d\n",($1+1)/2)}' `

MINEXP=$3
MAXEXP=$4
HALFEXP=$(echo "$MINEXP $MAXEXP" | awk '{printf("%d\n",$1+($2-$1)/2)}')

HEATMAPPER=$5

#First make the birdseye view of the figures
#export MAXCID
#export HALFCID # unnecessary to export...

# colormap string
# cluster IDs, then expression, then size of cluster.
sizemap="0:(255,255,255);200:(100,0,100)"
modmap="-2:(250,250,250);-1:(255,255,255);1:(77,172,38);$HALFCID:(245,245,245);$MAXCID:(208,28,139) 0:(0,0,200);$HALFEXP:(255,255,255);$MAXEXP:(200,0,0) $sizemap"




cat $OUTDIR/ordered_clusterset_means.txt| $HEATMAPPER -vC="$modmap" -vStrokeC="-" -vStrokeSC="black" -vFontSize=8 -vL="moduleIDs expr size"> $FIGDIR/ordered_clusterset_means.svg
#Now make the figures for each of the clusters of reasonable size
for FNAME in `ls $OUTDIR/clusterset*.txt`
do
	#echo $FNAME
	FIGNAME=${FNAME/.txt/.svg}
	FIGNAME=${FIGNAME/$OUTDIR}
	FIGNAME=$FIGDIR/$FIGNAME 	
 	cat $FNAME | $HEATMAPPER -vC="$modmap" -vStrokeC="-" -vStrokeSC="black" -vFontSize=8 -vL="moduleID expr size" > $FIGNAME
done

# Uncomment to make the all_assign.svg heatmap
#cat $OUTDIR/all_assign.txt | $HEATMAPPER -vC="$modmap" -vStrokeC="-" -vStrokeSC="black" -vFontSize=6 -vL="moduleID expr size" > $FIGDIR/all_assign.svg
