#! /usr/bin/env Rscript
# coding=utf-8

library(sys)
if (!require('devtools')) install.packages('devtools')
if (!require('PoreMeth2')) {
  options(timeout=9999999)
  devtools::install_github("Lab-CoMBINE/PoreMeth2")
}

library(PoreMeth2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
print("######################################################################################################################")
print("PoreMeth2 Arguments:")
print("Test:")
print(args[1])
print("Control:")
print(args[2])
print("Test name:")
print(args[3])
print("Control name:")
print(args[4])
print("Output:")
print(args[5])
print("omega:")
print(args[6])
print("eta")
print(args[7])
print("FW")
print(args[8])
print("AnnotationType")
print(args[9])
print("Annotate_Assembly")
print(args[10])
print("Statistics_Assembly")
print(args[11])
print("BetaThr")
print(args[12])
print("EntropyThr")
print(args[13])
print("PValueThr")
print(args[14])
print("AnalysisClass")
print(args[15])
print("betacovThr")
print(args[16])
print("######################################################################################################################")
print("reading files")

TableTest <- fread(args[1])
TableControl <- fread(args[2])

testname=args[3]
controlname=args[4]
output=args[5]
omeg=as.double(args[6])
et=as.double(args[7])
fw=as.integer(args[8])
AnnotType=args[9]
Annotate_Assembly=args[10]
Statistics_Assembly=args[11]
BetThr=as.double(args[12])
EntThr=as.double(args[13])
PThr=as.double(args[14])
AnalysisClass=args[15]
betacovThr=as.integer(args[16])

TableDMR_results=paste(output,"_DMRTable.txt",sep = "")
AnnotatedTableDMR_results=paste(output,"_DMRAnnotatedTable.txt",sep = "")
# PoreMeth2DMRStatistics_results=paste(output,"_DMRStatistics.txt",sep = "")
TestExpQualityPlot=paste(testname,"_ExpQualityPlot.png",sep = "")
ControlExpQualityPlot=paste(controlname,"_ExpQualityPlot.png",sep = "")
PairedExpQualityPlot=paste(output,"_paired_ExpQualityPlot.png",sep = "")

print("DMA process started")

# filt input data based on beta_cov
if (betacovThr != 0){
  TableTest <- TableTest[V6 >= betacovThr]
  TableControl <- TableControl[V6 >= betacovThr]
}

print(TableTest)
print(TableControl)

TableDMR  <-  PoreMeth2DMR(TableTest,TableControl, omega = omeg, eta = et, FW = fw)
write.table(TableDMR,TableDMR_results, sep="\t", row.names=F, quote=F)

# recommend to filter out DMRs with |$\Delta\beta$| < 0.2 to exclude unreliable results
print(class(TableDMR$DeltaBeta))
TableDMR$DeltaBeta <- as.numeric(TableDMR$DeltaBeta)
print(sum(is.na(TableDMR$DeltaBeta)))
TableDMR_filtered <- TableDMR[abs(TableDMR$DeltaBeta) < 0.2 ]

AnnotatedTableDMR <- PoreMethAnnotate2(TableDMR, NumProc = 5, AnnotationType = AnnotType, Assembly = Annotate_Assembly)
write.table(AnnotatedTableDMR,AnnotatedTableDMR_results, sep="\t", row.names=F, quote=F)

# PoreMeth2DMRStatistics(TableDMR, Assembly = Statistics_Assembly, BetaThr = BetThr, EntropyThr = EntThr, PValueThr = PThr)
# write.table(PoreMeth2DMRStatistics,PoreMeth2DMRStatistics_results, sep="\t", row.names=F, quote=F)

png(TestExpQualityPlot )
PoreMeth2SingleExpQualityPlot(TableTest)
dev.off()

png(ControlExpQualityPlot)
PoreMeth2SingleExpQualityPlot(TableControl)
dev.off()

png(PairedExpQualityPlot)
PoreMeth2PairedExpQualityPlot(TableTest,TableControl)
dev.off()
