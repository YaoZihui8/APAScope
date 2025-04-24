apalyzer <- function(out_dir) {
  library(GenomicRanges)
  library(GenomicFeatures)
  library("rtracklayer")
  library(dplyr)

  .annotatePASRegion <- function(pasdb, dbreigon, dfTXGeneMapping, flag)
  {
    x = unlist(dbreigon)
    x$IndexID = 1:length(x)
    x$transcript_id = names(x)
    names(x) = x$IndexID
    acc2sym <- data.frame(names(x), x$transcript_id)
    colnames(acc2sym) = c('IndexID', 'transcript_id')
    acc2sym = merge(acc2sym, dfTXGeneMapping, on = 'transcript_id', all.x = TRUE)
    x = x[acc2sym$IndexID]
    names(x) = acc2sym$gene_id
    dbreigon = split(x, names(x))
    dbreigon = GenomicRanges::reduce(dbreigon)
    olp = findOverlaps(pasdb, dbreigon)
    if (length(olp) > 0) {
      hit = pasdb[queryHits(olp),]
      hit$LOCATION = flag
      hit$PASid = paste(seqnames(hit), start(hit), strand(hit), sep = ":")
      hit = as.data.frame(hit)

    } else {
      hit = "NoHit"
    }
    return(hit)
  }

  .annotatePAS_V2 <- function(pasdb, TXDB, dfTXGeneMapping)
  {
    ##a.5UTR
    fiveUTRs = fiveUTRsByTranscript(TXDB, use.names = TRUE)
    fiveUTRhit = .annotatePASRegion(pasdb, fiveUTRs, dfTXGeneMapping, 'fiveUTR')
    ##b.intron
    introns = intronsByTranscript(TXDB, use.names = TRUE)
    intronhit = .annotatePASRegion(pasdb, introns, dfTXGeneMapping, 'intron')
    ##c.3UTR
    threeUTRs = threeUTRsByTranscript(TXDB, use.names = TRUE)
    threeUTRhit = .annotatePASRegion(pasdb, threeUTRs, dfTXGeneMapping, 'threeUTR')
    ##d.CDS
    cds = cdsBy(TXDB, by = "tx", use.names = TRUE)
    cdshit = .annotatePASRegion(pasdb, cds, dfTXGeneMapping, 'coding')
    ##e.all Exons
    exon = exonsBy(TXDB, by = "tx", use.names = TRUE)
    exonhit = .annotatePASRegion(pasdb, exon, dfTXGeneMapping, 'Exon')

    ## Creat header
    dbhit = as.data.frame(pasdb[1])
    dbhit$LOCATION = "NONE"
    dbhit$PASid = "NONE"
    for (hit in list(fiveUTRhit, intronhit, threeUTRhit, cdshit, exonhit)) {
      if (class(hit) == "data.frame") {
        dbhit = rbind(dbhit, hit)
      }
    }
    dbhit = dbhit[dbhit$PASid != "NONE",]
    finaldf = dbhit[, c('PASid', 'LOCATION', 'gene_id')]
    colnames(finaldf)[3] = "GENEID"
    finaldf$uniID1 = paste(finaldf$PASid, finaldf$GENEID, finaldf$LOCATION, sep = ":")
    finaldf$uniID2 = paste(finaldf$PASid, finaldf$GENEID, sep = ":")
    finaldf = finaldf[!duplicated(finaldf),]
    finaldf$ORDERID = 10
    if (nrow(finaldf[which(finaldf$LOCATION == 'fiveUTR'),]) > 0) {
      finaldf[which(finaldf$LOCATION == 'fiveUTR'),]$ORDERID = 1
    }
    if (nrow(finaldf[which(finaldf$LOCATION == 'intron'),]) > 0) {
      finaldf[which(finaldf$LOCATION == 'intron'),]$ORDERID = 2
    }
    if (nrow(finaldf[which(finaldf$LOCATION == 'coding'),]) > 0) {
      finaldf[which(finaldf$LOCATION == 'coding'),]$ORDERID = 3
    }
    if (nrow(finaldf[which(finaldf$LOCATION == 'threeUTR'),]) > 0) {
      finaldf[which(finaldf$LOCATION == 'threeUTR'),]$ORDERID = 4
    }
    finaldf = finaldf[with(finaldf, order(uniID2, ORDERID)),]
    finaldf = finaldf[!duplicated(finaldf$uniID2),]
    return(finaldf)
  }

  pas <- read.table(paste0(out_dir, "/dfPAS.txt"), header = FALSE)
  colnames(pas) <- c("seqnames", "start", "end", "width", "strand", "transcript_id", "gene_id", "gene_name", "PAS")
  dfPAS = as.data.frame(pas)
  dfTXGeneMapping = dfPAS[, c('transcript_id', 'gene_id')]
  dfTXGeneMapping = dfTXGeneMapping[!duplicated(dfTXGeneMapping),]
  dfTXGeneMapping = dfTXGeneMapping[!is.na(dfTXGeneMapping$transcript_id),]
  pasdb <- makeGRangesFromDataFrame(dfPAS, keep.extra.columns = TRUE)
  finaldf = .annotatePAS_V2(pasdb, TXDB, dfTXGeneMapping) #见下面的函数

  dfPAS$Chr = dfPAS$seqnames
  dfPAS$PASid = paste(dfPAS$seqnames, dfPAS$PAS, dfPAS$strand, sep = ":")
  dfinfo = dfPAS[, c('PASid', 'Chr', 'strand', 'PAS', 'gene_id', 'gene_name')]
  colnames(dfinfo) = c('PASid', 'Chr', 'strand', 'PAS', 'GENEID', 'gene_name')
  finaldf = merge(finaldf, dfinfo, by = c('PASid', 'GENEID'))

  # Extracting and filtering 3'UTR PASs
  library(dplyr)
  refUTRraw = .GTF2refUTRraw(TXDB, finaldf)

  .GTF2refUTRraw <- function(TXDB, finaldf)
  {
    ####### CDS end ############
    cds <- cdsBy(TXDB, "gene")
    x = unlist(cds)
    cds1 = split(x, names(x))
    test = as.data.frame(cds1)
    index = which(test$end - test$start > 1)
    allcds = test[index,]
    indexP = which(allcds$strand == '+')
    indexN = which(allcds$strand == '-')
    allcds$cdsend = ""
    allcds[indexP,]$cdsend = allcds[indexP,]$end
    allcds[indexN,]$cdsend = allcds[indexN,]$start
    MAX = aggregate(cdsend ~ group_name, allcds, max)
    MIN = aggregate(cdsend ~ group_name, allcds, min)
    names(MAX)[2] = 'max'
    names(MIN)[2] = 'min'
    maxmin = merge(x = MAX, y = MIN, by = "group_name")
    allcdssub = allcds[, c('group_name', 'strand')]
    allcdssub = allcdssub[!duplicated(allcdssub$group_name),]
    allcdsminmax = merge(x = allcdssub, y = maxmin, by = "group_name")
    indexP = which(allcdsminmax$strand == '+')
    indexN = which(allcdsminmax$strand == '-')
    allcdsminmax$cdsend = ""
    allcdsminmax[indexP,]$cdsend = allcdsminmax[indexP,]$max
    allcdsminmax[indexN,]$cdsend = allcdsminmax[indexN,]$min
    names(allcdsminmax)[1] = "GENEID"
    names(allcdsminmax)[2] = "Strand"
    allcdsminmax = allcdsminmax[, c('GENEID', 'Strand', 'cdsend')]


    ##### 3UTR PAS ######
    df3UTR = finaldf[which(finaldf$LOCATION == 'threeUTR'),]
    df3UTR = merge(df3UTR, allcdsminmax[, c('GENEID', 'cdsend')], by = c('GENEID'))
    ## Filter 3UTR PAS by CDS
    df3UTR = df3UTR[which(((df3UTR$strand == '+') & (df3UTR$PAS > df3UTR$cdsend)) | ((df3UTR$strand == '-') & (df3UTR$PAS < df3UTR$cdsend))),]

    df3UTR_P = df3UTR[which(df3UTR$strand == '+'),]
    df3UTR_P = df3UTR_P[with(df3UTR_P, order(GENEID, PAS)),]
    df3UTR_N = df3UTR[which(df3UTR$strand == '-'),]
    df3UTR_N = df3UTR_N[with(df3UTR_N, order(GENEID, -PAS)),]
    df3UTR = rbind(df3UTR_P, df3UTR_N)
    df3UTR_count = df3UTR %>%
      group_by(GENEID) %>%
      mutate(col = 1:n()) %>% ##create featureID
      mutate(count = n())
    df3UTR_count = data.frame(df3UTR_count)
    df3UTR_count$pA_type = 'M'
    if (nrow(df3UTR_count[which(df3UTR_count$col == 1 & df3UTR_count$count == 1),]) > 0) {
      df3UTR_count[which(df3UTR_count$col == 1 & df3UTR_count$count == 1),]$pA_type = 'S'
    }

    if (nrow(df3UTR_count[which(df3UTR_count$col == df3UTR_count$count & df3UTR_count$count > 1),]) > 0) {
      df3UTR_count[which(df3UTR_count$col == df3UTR_count$count & df3UTR_count$count > 1),]$pA_type = 'L'
    }

    if (nrow(df3UTR_count[which(df3UTR_count$col == 1 & df3UTR_count$count > 1),]) > 0) {
      df3UTR_count[which(df3UTR_count$col == 1 & df3UTR_count$count > 1),]$pA_type = 'F'
    }


    dfF = df3UTR_count[df3UTR_count$pA_type == 'F',][, c('GENEID', 'Chr', 'strand', 'PAS', 'gene_name')]
    dfL = df3UTR_count[df3UTR_count$pA_type == 'L',][, c('GENEID', 'Chr', 'strand', 'PAS', 'gene_name')]
    colnames(dfF) = c('GENEID', 'Chrom', 'Strand', 'First', 'gene_symbol')
    colnames(dfL) = c('GENEID', 'Chrom', 'Strand', 'Last', 'gene_symbol')
    dfFL = merge(dfF, dfL, by = c('GENEID', 'Chrom', 'Strand', 'gene_symbol'))
    dfFL = dfFL[!duplicated(dfFL$GENEID),]

    dfFLCDE = merge(dfFL, allcdsminmax, by = c('GENEID', 'Strand'))
    dfFLCDE = dfFLCDE[which(((dfFLCDE$Strand == '+') & (dfFLCDE$First > dfFLCDE$cdsend)) | ((dfFLCDE$Strand == '-') & (dfFLCDE$First < dfFLCDE$cdsend))),]
    dfFLCDE$First = as.integer(dfFLCDE$First)
    dfFLCDE$cdsend = as.integer(dfFLCDE$cdsend)
    dfFLCDE$distance = dfFLCDE$First - dfFLCDE$cdsend
    dfFLCDE[which(dfFLCDE$Strand == '-'),]$distance = dfFLCDE[which(dfFLCDE$Strand == '-'),]$cdsend - dfFLCDE[which(dfFLCDE$Strand == '-'),]$First

    index1 = which(dfFLCDE$Strand == '+' & dfFLCDE$distance < 100)
    dfFLCDE[index1,]$cdsend = dfFLCDE[index1,]$cdsend - (100 - dfFLCDE[index1,]$distance)
    index2 = which(dfFLCDE$Strand == '-' & dfFLCDE$distance < 100)
    dfFLCDE[index2,]$cdsend = dfFLCDE[index2,]$cdsend + (100 - dfFLCDE[index2,]$distance)
    refUTRraw = dfFLCDE[, c('gene_symbol', 'Chrom', 'Strand', 'First', 'Last', 'cdsend')]
    refUTRraw = refUTRraw[!is.na(refUTRraw$gene_symbol),]
    return(refUTRraw)
  }


  # Extracting IPAs
  dfIPAALL = .GTF2IPA(EDB, TXDB, finaldf)

  .GTF2IPA <- function(EDB, TXDB, finaldf)
  {
    ############ intron and SS ################
    introns = intronsByTranscript(TXDB, use.names = TRUE)
    dfintrons = as.data.frame(introns)
    colnames(dfintrons) = c('group', 'tx_id', 'Chrom', 'Start', 'End',
                            'width', 'Strand')
    # tx2gene <- mcols(transcripts(EDB, columns=c("tx_id", "gene_id",'gene_name')))
    # colnames(tx2gene)=c("tx_id", "GENEID",'gene_symbol')
    tx2gene = as.data.frame(EDB[which(EDB$type == 'transcript'),])[, c('transcript_id', 'gene_id', 'gene_name')]
    colnames(tx2gene) = c("tx_id", "GENEID", 'gene_symbol')
    tx2gene = tx2gene[!duplicated(tx2gene),]
    dfintrons = merge(dfintrons, tx2gene, by = 'tx_id')
    dfintrons = dfintrons[, c('GENEID', 'gene_symbol', 'Chrom', 'Start', 'End', 'Strand', 'width')]
    dfintrons = dfintrons[!duplicated(dfintrons),]

    ##ANNO 5SS and 3SS##
    dfintrons$SS5 = dfintrons$Start
    dfintrons$SS3 = dfintrons$End
    dfintrons[which(dfintrons$Strand == '-'),]$SS5 = dfintrons[which(dfintrons$Strand == '-'),]$End
    dfintrons[which(dfintrons$Strand == '-'),]$SS3 = dfintrons[which(dfintrons$Strand == '-'),]$Start

    ##Combine 5SS and 3SS##
    dfSS5 = dfintrons[, c("GENEID", "gene_symbol", "Chrom", "Strand", "SS5")]
    dfSS3 = dfintrons[, c("GENEID", "gene_symbol", "Chrom", "Strand", "SS3")]
    comname = c("GENEID", "gene_symbol", "Chrom", "Strand", "upstreamSS")
    colnames(dfSS5) = comname
    colnames(dfSS3) = comname
    dfSS5$type = 'SS5'
    dfSS3$type = 'SS3'
    dfupSS = rbind(dfSS5, dfSS3)
    dfupSS$SSID = paste0(dfupSS$GENEID, dfupSS$gene_symbol, dfupSS$Chrom, dfupSS$Strand, dfupSS$upstreamSS)
    dfupSS = dfupSS[!duplicated(dfupSS$SSID),]
    ############ build condidate downstream SS ##################
    dfdnSS3 = dfSS3
    dfdnSS3$SSID = paste0(dfdnSS3$GENEID, dfdnSS3$gene_symbol, dfdnSS3$Chrom, dfdnSS3$Strand, dfdnSS3$upstreamSS)
    dfdnSS3 = dfdnSS3[!duplicated(dfdnSS3$SSID),]
    colnames(dfdnSS3) = c("GENEID", "gene_symbol", "Chrom", "Strand", "downstreamSS", "type", "SSID")
    dfdnSS3$Chrom = paste0('chr', dfdnSS3$Chrom)
    dfupSS$Chrom = paste0('chr', dfupSS$Chrom)
    dfupSS = dfupSS[!is.na(dfupSS$gene_symbol),]
    dfdnSS3 = dfdnSS3[!is.na(dfdnSS3$gene_symbol),]

    ######### IPA and SS ############
    dfIPA = finaldf[(finaldf$LOCATION == 'intron'),]
    dfIPA = dfIPA[, c('PASid', 'Chr', 'PAS', 'strand', 'GENEID', 'gene_name')]
    colnames(dfIPA) = c('PASid', 'Chrom', 'Pos', 'Strand', 'GENEID', 'gene_symbol')
    dfIPA = dfIPA[!duplicated(dfIPA),]
    dfIPA = dfIPA[!is.na(dfIPA$gene_symbol),]
    DFIPASS = list(dfIPA, dfupSS, dfdnSS3)
    names(DFIPASS) = c('dfIPA', 'dfupSS', 'dfdnSS3')
    return(DFIPASS)
  }

  # Extracting 3' last exons
  dfIPAXXX = dfIPAALL$dfIPA
  dfupSSXXX = dfIPAALL$dfupSS
  dfdnSS3XXX = dfIPAALL$dfdnSS3

  old_text <- "chrchr"
  new_text <- "chr"
  dfupSSXXX_new <- apply(dfupSSXXX, 2, function(x) gsub(old_text, new_text, x))
  dfdnSS3XXX_new <- apply(dfdnSS3XXX, 2, function(x) gsub(old_text, new_text, x))
  dfupSSXXX <- dfupSSXXX_new
  dfdnSS3XXX <- dfdnSS3XXX_new
  dfupSSXXX <- as.data.frame(dfupSSXXX)
  dfupSSXXX$upstreamSS <- as.numeric(dfupSSXXX$upstreamSS)
  dfdnSS3XXX <- as.data.frame(dfdnSS3XXX)
  dfdnSS3XXX$downstreamSS <- as.numeric(dfdnSS3XXX$downstreamSS)

  dfLE_dfIPA = .GTF2LE(TXDB, dfIPAXXX, dfupSSXXX, dfdnSS3XXX)

  .GTF2LE <- function(TXDB, dfIPAXXX, dfupSSXXX, dfdnSS3XXX)
  {
    dfIPAXXX = as.data.frame(dfIPAXXX)
    dfupSSXXX = as.data.frame(dfupSSXXX)
    dfdnSS3XXX = as.data.frame(dfdnSS3XXX)
    ############ Hit upstream SS ##################
    dfHIT_up = merge(dfIPAXXX, dfupSSXXX, by = c("GENEID", "gene_symbol", "Chrom", "Strand"))
    dfHIT_up$dis = dfHIT_up$Pos - dfHIT_up$upstreamSS
    dfHIT_up[which(dfHIT_up$Strand == '-'),]$dis = dfHIT_up[which(dfHIT_up$Strand == '-'),]$upstreamSS - dfHIT_up[which(dfHIT_up$Strand == '-'),]$Pos
    dfHIT_up$HITID = paste0(dfHIT_up$GENEID, dfHIT_up$gene_symbol, dfHIT_up$PASid)
    dfHIT_up = dfHIT_up[which(dfHIT_up$dis > 100),]
    dfHIT_up = dfHIT_up[with(dfHIT_up, order(HITID, dis)),]
    dfHIT_up = dfHIT_up[!duplicated(dfHIT_up$HITID),]

    dfHIT_up = dfHIT_up[, c("GENEID", "gene_symbol", "Chrom", "Strand", "PASid", "Pos", "upstreamSS", "type", "dis")]
    colnames(dfHIT_up) = c("GENEID", "gene_symbol", "Chrom", "Strand", "PASid", "Pos", "upstreamSS", "type", "upstream_dis")

    ############ Hit downstream SS ##################
    dfHIT_dn = merge(dfIPAXXX, dfdnSS3XXX, by = c("GENEID", "gene_symbol", "Chrom", "Strand"))
    dfHIT_dn$dis = dfHIT_dn$downstreamSS - dfHIT_dn$Pos
    dfHIT_dn[which(dfHIT_dn$Strand == '-'),]$dis = dfHIT_dn[which(dfHIT_dn$Strand == '-'),]$Pos - dfHIT_dn[which(dfHIT_dn$Strand == '-'),]$downstreamSS
    dfHIT_dn$HITID = paste0(dfHIT_dn$GENEID, dfHIT_dn$gene_symbol, dfHIT_dn$PASid)
    dfHIT_dn = dfHIT_dn[which(dfHIT_dn$dis > 100),]
    dfHIT_dn = dfHIT_dn[with(dfHIT_dn, order(HITID, dis)),]
    dfHIT_dn = dfHIT_dn[!duplicated(dfHIT_dn$HITID),]

    dfHIT_dn = dfHIT_dn[, c("GENEID", "gene_symbol", "Chrom", "Strand", "PASid", "Pos", "downstreamSS", "dis")]
    colnames(dfHIT_dn) = c("GENEID", "gene_symbol", "Chrom", "Strand", "PASid", "Pos", "downstreamSS", "downstream_dis")

    ############ combine IPA SS ##################
    dfHIT_combine = merge(dfHIT_up, dfHIT_dn, by = c("GENEID", "gene_symbol", "Chrom", "Strand", "PASid", "Pos"))
    dfHIT_combine = dfHIT_combine[!duplicated(dfHIT_combine),]
    dfHIT_combine$IPA_up_start = dfHIT_combine$upstreamSS
    dfHIT_combine$IPA_up_end = dfHIT_combine$Pos
    dfHIT_combine$IPA_dn_start = dfHIT_combine$Pos
    dfHIT_combine$IPA_dn_end = dfHIT_combine$downstreamSS

    dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$IPA_up_start = dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$Pos
    dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$IPA_up_end = dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$upstreamSS
    dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$IPA_dn_start = dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$downstreamSS
    dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$IPA_dn_end = dfHIT_combine[which(dfHIT_combine$Strand == '-'),]$Pos

    dfIPAREF = dfHIT_combine[, c('PASid', "GENEID", 'gene_symbol', 'Chrom', 'Strand', 'Pos', 'upstreamSS', 'downstreamSS')]
    dfIPAREF = as.data.frame(dfIPAREF)
    dfIPAREF = dfIPAREF[!duplicated(dfIPAREF),]

    dfIPAsim = dfHIT_combine[, c('PASid', 'gene_symbol', 'Chrom', 'Strand', 'Pos', 'upstreamSS', 'downstreamSS')]
    dfIPAsim = as.data.frame(dfIPAsim)
    dfIPAsim = dfIPAsim[!duplicated(dfIPAsim),]
    dfIPAsim = dfIPAsim[!is.na(dfIPAsim$gene_symbol),]

    ############ TES ##############
    dfgenic = as.data.frame(genes(TXDB))
    dfgenic = dfgenic[, c(6, 1:3, 5)]
    names(dfgenic) <- c('GENEID', 'Chr', 'Start', 'End', 'Strand')
    #find the TES
    dfMinMAX = dfgenic
    dfMinMAX$TES = ''
    indexN = which(dfMinMAX$Strand == '-')
    indexP = which(dfMinMAX$Strand == '+')
    dfMinMAX[indexP,]$TES = dfMinMAX[indexP,]$End
    dfMinMAX[indexN,]$TES = dfMinMAX[indexN,]$Start
    dfTES = dfMinMAX[, c('GENEID', 'Chr', 'TES', 'Strand')]
    dfTES = dfTES[!duplicated(dfTES),]

    ###I exon###
    exondb = exonsBy(TXDB, by = "gene")
    exondb = GenomicRanges::reduce(exondb)
    exonsdf = as.data.frame(exondb)
    exonsdfP = exonsdf[which(exonsdf$strand == '+'),]
    exonsdfN = exonsdf[which(exonsdf$strand == '-'),]
    exonsdfN = exonsdfN[with(exonsdfN, order(group_name, -start)),]
    exonsdf = rbind(exonsdfP, exonsdfN)
    exonsdf = exonsdf[, -c(1)]

    data_EXrename = exonsdf %>%
      group_by(group_name) %>%
      mutate(col = 1:n()) %>% ##create exonID
      mutate(count = n())
    exonsdf = data.frame(data_EXrename)
    names(exonsdf)[7:8] = c('EXid', 'EXcount')

    dfLexon = exonsdf[which(exonsdf$EXid == exonsdf$EXcount),]
    dfLexon$LEstart = dfLexon$start
    dfLexon[which(dfLexon$strand == '-'),]$LEstart = dfLexon[which(dfLexon$strand == '-'),]$end
    dfLexon = dfLexon[, c('group_name', 'seqnames', 'strand', 'LEstart')]
    colnames(dfLexon) = c('GENEID', 'Chr', 'Strand', 'LEstart')
    dfLE = merge(dfLexon, dfTES, by = c('GENEID', 'Chr', 'Strand'))
    dfLE = dfLE[!duplicated(dfLE),]

    dfIPAused = dfIPAREF[, c("GENEID", 'gene_symbol')]
    dfIPAused = dfIPAused[!duplicated(dfIPAused),]

    dfLE = merge(dfLE, dfIPAused, by = 'GENEID')
    dfLE = dfLE[!duplicated(dfLE),]
    dfLE$Chr = paste0('chr', dfLE$Chr)
    dfLE = dfLE[, c('gene_symbol', 'Chr', 'Strand', 'LEstart', 'TES')]
    colnames(dfLE) = c('gene_symbol', 'Chrom', 'Strand', 'LEstart', 'TES')
    dfLE = dfLE[!is.na(dfLE$gene_symbol),]
    dfLE_dfIPA = list(dfIPAsim, dfLE)
    names(dfLE_dfIPA) = c('dfIPAsim', 'dfLE')
    return(dfLE_dfIPA)
  }

  # Finalizing references
  PASREF = list(refUTRraw, dfLE_dfIPA$dfIPAsim, dfLE_dfIPA$dfLE)
  names(PASREF) = c('refUTRraw', 'dfIPA', 'dfLE')

  refUTRraw = PASREF$refUTRraw
  dfIPAraw = PASREF$dfIPA
  dfLEraw = PASREF$dfLE

  write.table(refUTRraw, file = paste0(out_dir,"UTRraw.tsv"), sep = "\t", row.names = FALSE)
  write.table(dfIPAraw, file = paste0(out_dir,"IPAraw.tsv"), sep = "\t", row.names = FALSE)
  write.table(dfLEraw, file = paste0(out_dir,"LEraw.tsv"), sep = "\t", row.names = FALSE)


}