echo "Making plots"
Rscript resources/genetics/plot_gwas.R \
	${section_12_dir}/cellcount_${batch}.loco.mlma.gz \
	10 \
	1 \
	3 \
 	2 \
	TRUE \
	0 \
	0 \
	0 \
	0 \

echo "Successfully performed GWAS for cell type ${batch}"
