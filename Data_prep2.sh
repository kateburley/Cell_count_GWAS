# Create covariates for cellcounts GWAS
# FID IID Smoking Age_numeric Sex

covnames <- c("Age_numeric", "Sex_factor")
covnames <- covnames[covnames %in% names(allcovs)]

dat <- merge(fam, smok, by.x="V2", by.y="IID")
dat <- subset(dat, select=c("V1", "V2", "Smoking"))
if(length(covnames) > 0)
{
	dat <- merge(dat, subset(allcovs, select=c("IID", covnames)), by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", "Smoking", covnames))
}
names(dat)[names(dat) == "Smoking"] <- "Smoking_numeric"
write_covs(dat, paste0(out_file, ".cellcounts"))
