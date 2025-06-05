echo "Making plots"
Rscript Cell_count_gwas/plot_gwas2.R \
	${section_12_dir}/cellcount_${batch}.loco.mlma.gz \
	9 \
	1 \
	3 \
 	2 \
	TRUE \
	0 \
	0 \
	0 \
	0 \

echo "Successfully performed GWAS for cell type ${batch}"
