
#Prepare Q-motif data from example

#map motifs to genes with matchMotifToGenePerTF2.py
python aggregateSignalMotifNet/matchMotifToGenePerTF2.py example_files/mouse_chr1_motif_example.txt example_files/Mus_musculus.NCBIM37.67.TSS.txt 2500 2500 example_files/mouse_chr1_example.txt

export exeMN=aggregateSignalMotifNet/aggregateSignalMotifNet

#aggregate signal for these mapped motif instances for each of the fiven example data sets
for S in ESC MEF MEF_48hr pips1 pips2
do
    eval ${exeMN} example_files/mouse_chr1_example.txt example_files/mm9.fa.fai example_files/mouse_${S}_chr1_example.counts example_files/mouse_${S}_chr1_example_Q_motif -n 1
done

#prepare list of aggregated value data files
ls example_files/*Q_motif_aggregated_values.txt > example_files/mouse_chr1_example_list.txt

#merge the aggregated motif values.
eval mergedata/mergeData example_files/mouse_chr1_example_list.txt example_files/mouse_chr1_example_merged.txt

#remove rows that are not uniform, that are associated with missing or no data in any data set.
grep -F -v nodata example_files/mouse_chr1_example_merged.txt > example_files/mouse_chr1_example_merged_uniform.txt

#apply log transformation and quantile normalization
matlab -nodisplay -r "quantile_normalize('example_files/mouse_chr1_example_merged_uniform.txt','example_files/mouse_chr1_example_merged_normalized.txt','doLog');exit;"


#Prepare promoter signal data from example
exeR=aggregateSignalRegion_nonLog/aggregateSignal

#prepare input set of promoter regions
awk '$2=="chr1"{printf("%s\t%i\t%i\t%s\t+\n",$2,$3-2500,$4+2500,$1)}' example_files/Mus_musculus.NCBIM37.67.TSS.txt > example_files/promoters.bed

#aggregate signal for promoter regions from each condition/data set
for S in ESC MEF MEF_48hr pips1 pips2
do
    eval ${exeR} example_files/promoters.bed example_files/mm9.fa.fai example_files/mouse_${S}_chr1_example.counts example_files/mouse_${S}_chr1_example_promoter -n 1
done

#merge and normalize the promoter signal data as for the Q-motif signal data
#list of aggregated promoter signal files from aggregateSignalRegion
ls example_files/*promoter_aggregate.txt > example_files/mouse_chr1_example_promoter_list.txt
#merge promoter signal data across conditions
eval mergedata/mergeData example_files/mouse_chr1_example_promoter_list.txt example_files/mouse_chr1_example_promoter_merged.txt
#select data for regions uniformly represented in the data.
grep -F -v nodata example_files/mouse_chr1_example_promoter_merged.txt > example_files/mouse_chr1_example_promoter_merged_uniform.txt
#apply log transformation and quantile nromalization
matlab -nodisplay -r "quantile_normalize('example_files/mouse_chr1_example_promoter_merged_uniform.txt','example_files/mouse_chr1_example_promoter_merged_normalized.txt','doLog');exit;"

#Prepare feature input data files from example

#First get the number of regulatory features, so the number of motifs +1, for ATAC signal data
NR=`awk '{print $1}' example_files/mouse_chr1_example_merged_normalized.txt | sed 's/_/\t/g' | cut -f1 | sort -u | wc -l | awk '{print $1+1}'`
#Second get number of genes in either the Q-motif or promoter signal data.
NG=`cat <(awk '{print $1}' example_files/mouse_chr1_example_merged_normalized.txt | sed 's/_/\t/g' | cut -f2 | sort -u) <(cut -f1 example_files/mouse_chr1_example_promoter_merged_normalized.txt) | sort -u | wc -l`

c=2
for S in ESC MEF MEF_48hr pips1 pips2
do
    export S=${S}
    #add header line with number of features and genes
    printf "${NR}\t${NG}\n" > example_files/mouse_${S}_chr1_example_features.txt
    #add Q-motif data 
    awk -v c=$c '{split($1,a,"_");printf("%s\t%s_%s\t%s\n",a[1],a[2],ENVIRON["S"],$c)}' example_files/mouse_chr1_example_merged_normalized.txt >> example_files/mouse_${S}_chr1_example_features.txt
    #add promoter signal data
     awk -v c=$c '{printf("ATAC\t%s_%s\t%s\n",$1,ENVIRON["S"],$c)}' example_files/mouse_chr1_example_promoter_merged_normalized.txt >> example_files/mouse_${S}_chr1_example_features.txt
    c=$((c+1))
done
