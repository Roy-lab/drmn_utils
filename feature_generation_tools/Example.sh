exeMN=aggregateSignalMotifNet/aggregateSignalMotifNet

#python aggregateSignalMotifNet/matchMotifToGenePerTF2.py example_files/mouse_chr1_motif_example.txt example_files/Mus_musculus.NCBIM37.67.TSS.txt 2500 2500 example_files/mouse_chr1_example.txt

#for S in ESC MEF MEF_48hr pips1 pips2
#do
#	eval ${exeMN} example_files/mouse_chr1_example.txt example_files/mm9.fa.fai example_files/mouse_${S}_chr1_example.counts example_files/mouse_${S}_chr1_example_Q_motif -n 1
#done

#ls example_files/*Q_motif_aggregated_values.txt > example_files/mouse_chr1_example_list.txt 

#eval mergedata/mergeData example_files/mouse_chr1_example_list.txt example_files/mouse_chr1_example_merged.txt

#grep -F -v nodata example_files/mouse_chr1_example_merged.txt > example_files/mouse_chr1_example_merged_uniform.txt

#matlab -nodisplay -r "quantile_normalize('example_files/mouse_chr1_example_merged_uniform.txt','example_files/mouse_chr1_example_merged_normalized.txt','doLog');exit;"

exeR=aggregateSignalRegion_nonLog/aggregateSignal

awk '$2=="chr1"{printf("%s\t%i\t%i\t%s\t+\n",$2,$3-2500,$4+2500,$1)}' example_files/Mus_musculus.NCBIM37.67.TSS.txt > example_files/promoters.bed

for S in ESC MEF MEF_48hr pips1 pips2
do
    eval ${exeR} example_files/promoters.bed example_files/mm9.fa.fai example_files/mouse_${S}_chr1_example.counts example_files/mouse_${S}_chr1_example_promoter -n 1
done

ls example_files/*promoter_aggregate.txt > example_files/mouse_chr1_example_promoter_list.txt

eval mergedata/mergeData example_files/mouse_chr1_example_promoter_list.txt example_files/mouse_chr1_example_promoter_merged.txt

grep -F -v nodata example_files/mouse_chr1_example_promoter_merged.txt > example_files/mouse_chr1_example_promoter_merged_uniform.txt

matlab -nodisplay -r "quantile_normalize('example_files/mouse_chr1_example_promoter_merged_uniform.txt','example_files/mouse_chr1_example_promoter_merged_normalized.txt','doLog');exit;"

NR=`awk '{print $1}' example_files/mouse_chr1_example_merged_normalized.txt | sed 's/_/\t/g' | cut -f1 | sort -u | wc -l | awk '{print $1+1}'`
NG=`awk '{print $1}' example_files/mouse_chr1_example_merged_normalized.txt | sed 's/_/\t/g' | cut -f2 | sort -u | wc -l`

c=2
for S in ESC MEF MEF_48hr pips1 pips2
do
	export S=${S}
    printf "${NR}\t${NG}\n" > example_files/mouse_${S}_chr1_example_features.txt
    awk -v c=$c '{split($1,a,"_");printf("%s\t%s_%s\t%s\n",a[1],a[2],ENVIRON["S"],$c)}' example_files/mouse_chr1_example_merged_normalized.txt >> example_files/mouse_${S}_chr1_example_features.txt
	 awk -v c=$c '{printf("ATAC\t%s_%s\t%s\n",$1,ENVIRON["S"],$c)}' example_files/mouse_chr1_example_promoter_merged_normalized.txt >> example_files/mouse_${S}_chr1_example_features.txt
    c=$((c+1))
done



