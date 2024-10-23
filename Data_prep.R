# Create covariates files for cellcounts GWAS
# FID IID Smoking Age_numeric Sex

write_covs <- function(dat, filename)
{
	fcols <- grep("_factor", names(dat))
	ncols <- grep("_numeric", names(dat))
	fdat <- data.frame(dat[,1:2], dat[,fcols])
	ndat <- data.frame(dat[,1:2], dat[,ncols])
	write.table(fdat, file=paste0(filename, ".factor"), row=F, col=F, qu=F)
	write.table(ndat, file=paste0(filename, ".numeric"), row=F, col=F, qu=F)
	write.table(dat, file=filename, row=F, col=F, qu=F)
}

arguments <- commandArgs(T)

covs_file <- arguments[1]
aar_file <- arguments[2]
smok_file <- arguments[3]
fam_file <- arguments[4]
out_file <- arguments[5]
covs_orig <- arguments[6]

allcovs <- read.table(covs_file, he=T, stringsAsFactors=FALSE)
smok <- read.table(smok_file, he=T, stringsAsFactors=FALSE)
fam <- read.table(fam_file, stringsAsFactors=FALSE)[,1:2]

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



# Create cell counts files for cellcounts GWAS and do inverse normal transformation
