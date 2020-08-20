## cd ../results
#nohup bash ../scripts/snakebash.sh $spec > nohup_${spec}.txt &
#augur refine  --tree tree_pruned_2020-06-25_b1b2Na.nwk \
#              --alignment .fasta             --metadata ../data/metadata_2020-06-25b1b2N.tsv    --output-tree tree_20
#20-06-25_b1b2N.nwk             --output-node-data branch_lengths_2020-06-25_b1b2N.json             --timetre
#e             --coalescent opt


#spec=2020-06-25_b1-3QCNR_countryMax300
#spec=2020-06-25_b1-3QCNR_countryMax300_1k
#spec=2020-06-25_b1-3QCNR_countryMax300_k100
spec=$1

#spec=Gisaid0625_b1b2_k1k_countryMax300
#spec=2020-06-25_b1b2N_28K
tree0=tree_raw_${spec}.nwk
tree=tree_${spec}.nwk
ali=aligned_${spec}.fasta
branch=branch_lengths_${spec}.json
nt=nt_muts_${spec}.json
aa=aa_muts_${spec}.json
traits=traits_${spec}.json
ref="../config/reference.gb"
metadata="../data/metadata_2020-06-25_b1-3QC.tsv"
#metadata="../data/metadata_${spec}.tsv"
branches="branch_lengths_${spec}.json"
colors="../config/color_schemes.tsv"
auspice_config="../config/auspice_config.json"
lat_longs="../config/lat_long.tsv"
output="../auspice/covid19uae_${spec}.json"

augur tree --alignment aligned_${spec}.fasta --output tree_raw_${spec}.nwk --nthreads 60 

augur refine  --tree $tree0 --alignment $ali --output-tree $tree --metadata $metadata --output-node-data $branch --timetree --coalescent opt --root MN908947
## ancestral
augur ancestral --tree $tree --alignment $ali --output-node-data $nt --inference "joint"

augur translate --tree $tree --reference-sequence $ref --ancestral-sequences $nt --output-node-data $aa

augur traits --tree $tree --metadata $metadata --output-node-data $traits --columns "country" --confidence

augur export v2 --tree $tree --metadata $metadata --node-data $branches $traits $nt $aa --colors $colors --lat-longs $lat_longs --auspice-config $auspice_config --output $output

 
