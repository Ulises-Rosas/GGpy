# for aln in *.r2_para_no_lavaretus_TBL_tlike_aln;do
#     new_name=$( echo $aln | sed "s/r2_para_no_lavaretus_TBL_tlike_aln/fasta/g" );
#     mv $aln $new_name;
# done


for tree in *.r2_para_no_lavaretus_TBL_tlike_aln.nex.treefile;do
    new_name=$( echo $tree | sed "s/r2_para_no_lavaretus_TBL_tlike_aln.nex.treefile/tree/g" );
    mv $tree $new_name;
done
