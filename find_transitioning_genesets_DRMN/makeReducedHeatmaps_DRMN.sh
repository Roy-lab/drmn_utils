#!/bin/bash
# Plots a manually reduced bird's eye view
set -u
if [ $# -lt 4 ]
then
	echo "Too few arguments $#"
	echo "Usage: ./makeTransitionHeatmaps.sh myfile numberofclusters minExpr maxExpr"
	exit
fi
MYFILE=$(readlink -f $1)
mybase=$(basename $MYFILE)
OUTDIR=${MYFILE/$mybase}

#MAXCID=`expr $2 - 1` # actually, these are getting offset by 1
MAXCID=$2
HALFCID=`echo $MAXCID| awk '{printf("%d\n",($1+1)/2)}' `

MINEXP=$3
MAXEXP=$4
HALFEXP=$(echo "$MINEXP $MAXEXP" | awk '{printf("%d\n",$1+($2-$1)/2)}')

HEATMAPPER=/mnt/dv/wid/projects2/Roy-common/programs/scripts/figscripts/Heatmap.awk

# colormap string
# cluster IDs, then expression, then size of cluster.
sizemap="0:(255,255,255);200:(100,0,100)"
modmap="-2:(250,250,250);-1:(255,255,255);1:(77,172,38);$HALFCID:(245,245,245);$MAXCID:(208,28,139) 0:(0,0,200);$HALFEXP:(255,255,255);$MAXEXP:(200,0,0) $sizemap 2:(255,255,255);5:(0,100,100)"

cat $MYFILE| $HEATMAPPER -vC="$modmap" -vStrokeC="-" -vStrokeSC="black" -vFontSize=8 -vL="moduleIDs expr size -log10qval" -vD="-2"> $OUTDIR/${mybase}.svg

convert ${OUTDIR}/${mybase}.svg ${OUTDIR}/${mybase}.png
