###<Script provided by Pio, which was then manually checked and verified by me>

#!/usr/bin/env Rscript

options(warn = 0)

## Libraries

if (suppressPackageStartupMessages(!require("this.path", quietly = TRUE))) {
  install.packages("this.path")
}

if (suppressPackageStartupMessages(!require("data.table", quietly = TRUE))) {
  install.packages("data.table")
}

if (suppressPackageStartupMessages(!require("parallel", quietly = TRUE))) {
  install.packages("parallel")
}

if (suppressPackageStartupMessages(!require("LncFinder", quietly = TRUE))) {
  lx("Intalling package [LncFinder]")
  install.packages("LncFinder", verbose =F) #
}

if (suppressPackageStartupMessages(!require("qualV", quietly = TRUE))) {
  lx("Intalling package [qualV]")
  install.packages("qualV", verbose =F) #
}

if (suppressPackageStartupMessages(!require("stringi", quietly = TRUE))) {
  lx("Intalling package [stringi]")
  install.packages("stringi", verbose =F) #
}


# Reads a fasta file
ffasta <- function(f) {
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    return(data.table(name = as.character(), seq = as.character()))
  } else {
    fa_raw <- fread(f, header = F, fill = T, sep = "\n")
    fa_raw[, h := grepl(">", V1)]
    fa_oneline <- fa_raw[, .(paste0(V1, collapse = "")), by = rleid(h)]
    return(data.table(name = fa_oneline[rep(c(TRUE, FALSE), length = .N)]$V1, 
                      seq = fa_oneline[rep(c(FALSE, TRUE), length = .N)]$V1))
  }
}

## Writes a fasta file
wfasta <- function(fasta_data, f) {
  fileconn <- file(f)
  writeLines(t(fasta_data), fileconn)
  close(fileconn)
}

### Detection of small, degraded TIR
short_tir <- function(a,b) {
  lcs <- LCS(strsplit(a,"")[[1]], strsplit(b,"")[[1]])
  return(nchar(paste0(lcs$LCS,collapse = ""))>7)
}

rc <- function(s) {
  stri_reverse(chartr("ACGTacgt","TGCAtgca",s))
}

three_orfs <- function(seq) {
  orfs <- rbindlist(find_orfs(seq, reverse.strand = TRUE, 
                              max.only = F), fill = T)
  orfs_list <- orfs[,first(ORF.Len), ORF.Stop]
  orfs_list_ordered <- orfs_list[order(-V1)]$V1
  return(data.table(orf1 = orfs_list_ordered[1],
                    orf2 = orfs_list_ordered[2],
                    orf3 = orfs_list_ordered[3]))
}


f <- commandArgs(trailingOnly=TRUE)[1]
tes <- ffasta(f)
tes[,seq:=toupper(seq)]
tes[,lente:= nchar(seq),name]



temp <- paste0("temp",make.names(gsub(".*/","",f)))
dir.create(temp,showWarnings = F)
setwd(temp)


orfs_tes <- rbindlist(mclapply(tes$seq,three_orfs))
tes <- cbind(tes,orfs_tes)
wfasta(tes[, c("name","seq")], paste0("tmp-tes"))
system(paste0("makeblastdb -in tmp-tes -dbtype nucl 1> /dev/null"))
te_data <- fread(cmd= paste0("blastn -query tmp-tes -db tmp-tes -task blastn -num_threads ", 
                             detectCores(), 
                             " -evalue 500 -outfmt 6 -word_size 11 -gapopen 4 -gapextend 1 -reward 1 -penalty -1"), header = F)
if (nrow(te_data) > 0) {
  colnames(te_data) <-c("qseqid","sseqid", "pident" , "length","mismatch", 
                        "gapopen","qstart","qend","sstart","send", 
                        "evalue", "bitscore")
  te_data_mix <- te_data[qseqid != sseqid] # Non self matches, to deal with later
  te_data_rep <- te_data[qseqid == sseqid & length >20 & 
                           (qstart != sstart | qend != send)][,.(N=.N,maxRep=max(length)),qseqid] # large self matches that are not trivial. To filter with MaxRep.
  te_data  <- te_data[!(qseqid == sseqid & qstart == sstart & qend == send)] # Self matches that are not the whole element, to check for TIR and LTR
  ### TIR, LTR detection
  te_data <- te_data[qseqid == sseqid]
  te_data[,len:=abs(qend-qstart)]
  te_data <- te_data[order(qstart)][order(-len)]
  te_data <- merge(te_data, 
                   tes[,c("name","lente", "orf1", "orf2", "orf3")]
                   [,name:=gsub(">","",name)], by.x = "qseqid", by.y="name")
  te_data[,lgap:=.(min(qstart,qend,sstart,send)-1), by=1:nrow(te_data)]
  te_data[,rgap:=.(lente-max(qstart,qend,sstart,send)), by=1:nrow(te_data)]
  te_data[,check:=((qend+qstart-lente)/2) * ((send+sstart-lente)/2)<0, 
          by=1:nrow(te_data)]  ### Is each match at one side of the center?
  te_data <- te_data[check==T]
  te_data[,check_type:=((qend-qstart)*(send-sstart))<0, by=1:nrow(te_data)]
  te_data[,tgap :=lgap+rgap]
  te_data <- te_data[order(tgap)][,.SD[1],qseqid]
  te_data$type <- "LTR"
  te_data[check_type==T,type:= "TIR"]
  ### TIR, LTR detection END
  tes <- merge(tes[,name:=gsub(">","",name)], 
               te_data[,c("qseqid","length","lgap","rgap", "type")], 
               by.x = "name", by.y = "qseqid", all.x=T)
  tes[is.na(type), rgap :=0]
  tes[is.na(type), lgap :=0]
  tes[is.na(type), length :=0]
  
  tes <- merge(tes,te_data_rep, by.x = "name", by.y = "qseqid", all.x=T)
  tes_count <- nrow(tes)
  tes$pas <- unlist(lapply(stri_locate_all_regex(tes$seq,paste0(strrep("A",10),"|",strrep("T",10))),function(x) {min(unlist(x))-1}))
  tes$eas <- nchar(tes$seq)-unlist(lapply(stri_locate_all_regex(tes$seq,paste0(strrep("A",10),"|",strrep("T",10))),function(x) {max(unlist(x))}))

  tes[,pa:=min(pas,eas), by=.I]
  tes[substr(seq,1,5)=="AAAAA" | substr(seq,1,5)=="TTTTT" | substr(seq,lente-4,lente)=="AAAAA" | substr(seq,lente-4,lente)=="TTTTT" ,pa:=0, by = .I]

  tes[,st:=short_tir(substr(seq,1,10), rc(substr(seq,nchar(seq)-9,nchar(seq)))), .I]
  tes[st==T & is.na(type), type:="TIR"]
  tes[st==T & is.na(type), lgap := 0]
  tes[st==T & is.na(type), rgap := 0]
  tes[st==T & is.na(type), length:=8]
  
}
setwd("..")
unlink(temp, recursive = TRUE)
tes <- tes[order(-lente)]
stats_data <- tes[,c("name","lente","type",
                     "length","lgap","rgap","pa","orf1","orf2","orf3","maxRep")]


### Last fixes
stats_data$pa <- stats_data$pa<10
stats_data[is.na(pa),pa:=F]
colnames(stats_data) <- c("name","TE_len","struct_type","struct_len","left_gap",
                          "right_gap","polyA","orf1","orf2","orf3","maxTR")

fwrite(stats_data, paste0(f, ".stats"), 
       quote =  F, row.names = F, sep ="\t")
