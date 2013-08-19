GAPDH_1095_PMT0 <- rppa.load(baseUrl="http://10.149.64.8:8080/MIRACLE/spotExport/", barcode="1095")
TFR_1096_PMT0 <- rppa.load(baseUrl="http://10.149.64.8:8080/MIRACLE/spotExport/", barcode="1096")
CD44_1097_PMT0 <- rppa.load(baseUrl="http://10.149.64.8:8080/MIRACLE/spotExport/", barcode="1097")
CDH1_1098_PMT0 <- rppa.load(baseUrl="http://10.149.64.8:8080/MIRACLE/spotExport/", barcode="1098")
PTEN_1099_PMT0 <- rppa.load(baseUrl="http://10.149.64.8:8080/MIRACLE/spotExport/", barcode="1099")

MCF7.list <- list(GAPDH_1095_PMT0, TFR_1096_PMT0, CD44_1097_PMT0, CDH1_1098_PMT0, PTEN_1099_PMT0)
MCF7.list.sdc <- foreach(slide=MCF7.list)%do% { rppa.serialDilution(subset(slide, NumberOfCellsSeeded==2000), select.columns.A="Treatment", select.columns.B="NumberOfCellsSeeded", select.columns.fill="CellLine")}
mydata <- foreach(sdc=MCF7.list.sdc, .combine=cbind) %do% { 
  x <- log2(sdc$concentrations) 
  x <- x - median(x, na.rm=T)
}
colnames(mydata) <- c("GAPDH", "TFR", "CD44", "CDH1", "PTEN")

#
# Housekeeping
mydata.house <- sweep(mydata,1,mydata[,1],"-")

#
# Median
rowmedianB <- apply(mydata,1,median,na.rm=T)
colmedianB <- apply(mydata,2,median,na.rm=T)
temp <- sweep(mydata,2,colmedianB)
mydata.med <- sweep(temp,1,rowmedianB)
mydata.med <- sweep(mydata,1,rowmedianB)

#
# Variable Slope
colMedianB <- apply(mydata,2,median,na.rm=T)
rowMedianB <- apply(mydata,1,median,na.rm=T)
xij <- sweep(mydata,2,colMedianB)
gammaB <- estimateGamma(xij,method="other")
mydata.gam <- sweep(xij,2,gammaB,"/")
mydata.gam <- sweep(mydata.gam,1,rowMedianB,"-")

mydata.gam.melt <- melt(cbind(MCF7.list.sdc[[1]][,c("Sample", "A", "B")], mydata.gam), variable_name="Fill")
mydata.gam.melt$concentrations <- 2^mydata.gam.melt$value
mydata.gam.melt$upper <- 2^mydata.gam.melt$value
mydata.gam.melt$lower <- 2^mydata.gam.melt$value

rppa.proteinConc.plot(mydata.gam.melt, scales="free", error.bars=F, title="Variable slope normalization")