#!/bin/bash

set -e
source ./config

# Check cell type variable is valid
batch=${1}
re='^[0-9]+$'
ncellcounts=`wc -l ${home_directory}/results/03/cellcounts_summary.txt | awk '{ print $1 }'`
if ! [[ $batch =~ $re ]] ; then
	echo "error: Cell type variable is not valid"
	echo "Please provide a number between 1 and ${ncellcounts}"
	echo "Usage: ${0} [cell type]"
	exit 1
fi

if [ "${batch}" -gt "${ncellcounts}" ]; then
	echo "error: Cell type variable is not valid"
	echo "Please provide a number between 1 and ${ncellcounts}"
	echo "Usage: ${0} [cell type]"
	exit 1
fi

if [ "${batch}" -lt "1" ]; then
	echo "error: Cell type variable is not valid"
	echo "Please provide a number between 1 and ${ncellcounts}"
	echo "Usage: ${0} [cell type]"
	exit 1
fi

exec &> >(tee ${home_directory}/results/12/logs${batch})
print_version

# Generate GWAS covariates plus INT cell counts
echo "Generating GWA covariates"
Rscript resources/cellcounts/create_cellcount_GWAS_files.R \
	${covariates_combined}.txt \
	${smoking_pred}.txt \
	${bfile}.fam \
	${gwas_covariates} \
	${covariates}
 
age=`awk '{print $4}' <${gwas_covariates}.cellcounts.numeric |sort -u |wc -l`
echo "Age variable has $age levels"

sex=`awk '{print $3}' <${gwas_covariates}.cellcounts.factor |sort -u |wc -l`
echo "Sex variable has $sex levels"

n23=`grep ^23 ${bfile}.bim | wc -l`

if [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${home_directory}/processed_data/methylation_data/all_covariates.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--covar ${gwas_covariates}.cellcounts.factor \
	--autosome \
	--out ${home_directory}/results/12/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}
fi

if [ "$n23" -gt "0" ] && [ "$sex" -gt "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--chr 23 \
	--mlma \
	--pheno ${home_directory}/processed_data/methylation_data/all_covariates.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--covar ${gwas_covariates}.cellcounts.factor \
	--grm ${grmfile_all} \
	--out ${section_12_dir}/cellcount_${batch}_chr23 \
	--thread-num ${nthreads} \
	--mpheno ${batch}

fi


if [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${home_directory}/processed_data/methylation_data/all_covariates.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--autosome \
	--out ${section_12_dir}/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}
fi

if [ "$n23" -gt "0" ] && [ "$sex" -eq "1" ] && [ "$age" -gt "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--chr 23 \
	--mlma \
	--pheno ${home_directory}/processed_data/methylation_data/all_covariates.plink \
	--qcovar ${gwas_covariates}.cellcounts.numeric \
	--grm ${grmfile_all} \
	--out ${section_12_dir}/cellcount_${batch}_chr23 \
	--thread-num ${nthreads} \
	--mpheno ${batch}

	head -n1 ${section_12_dir}/cellcount_${batch}.loco.mlma >${section_12_dir}/cellcount_${batch}.loco
	tail -q -n +2 ${section_12_dir}/cellcount_${batch}.loco.mlma ${section_12_dir}/cellcount_${batch}_chr23.mlma >>${section_12_dir}/cellcount_${batch}.loco
	mv ${section_12_dir}/cellcount_${batch}.loco ${section_12_dir}/cellcount_${batch}.loco.mlma
	rm ${section_12_dir}/cellcount_${batch}_chr23.mlma

fi

if [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--mlma-loco \
	--pheno ${home_directory}/processed_data/methylation_data/all_covariates.plink \
	--covar ${gwas_covariates}.cellcounts.factor \
	--autosome \
	--out ${section_12_dir}/cellcount_${batch} \
	--thread-num ${nthreads} \
	--mpheno ${batch}
fi

if [ "$n23" -gt "0" ] && [ "$sex" -gt "1" ] && [ "$age" -eq "1" ]
then
${gcta} \
	--bfile ${bfile} \
	--chr 23 \
	--mlma \
	--pheno ${home_directory}/processed_data/methylation_data/all_covariates.plink \
	--covar ${gwas_covariates}.cellcounts.factor \
	--grm ${grmfile_all} \
	--out ${section_12_dir}/cellcount_${batch}_chr23 \
	--thread-num ${nthreads} \
	--mpheno ${batch}

	head -n1 ${section_12_dir}/cellcount_${batch}.loco.mlma >${section_12_dir}/cellcount_${batch}.loco
	tail -q -n +2 ${section_12_dir}/cellcount_${batch}.loco.mlma ${section_12_dir}/cellcount_${batch}_chr23.mlma >>${section_12_dir}/cellcount_${batch}.loco
	mv ${section_12_dir}/cellcount_${batch}.loco ${section_12_dir}/cellcount_${batch}.loco.mlma
	rm ${section_12_dir}/cellcount_${batch}_chr23.mlma

fi



echo "Compressing results"
gzip -f ${section_12_dir}/cellcount_${batch}.loco.mlma

echo "Making plots"
Rscript resources/genetics/plot_gwas.R \
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
