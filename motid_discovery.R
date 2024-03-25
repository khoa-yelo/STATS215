BiocManager::install("rGADEM")

library(rGADEM)
library(BSgenome.Hsapiens.UCSC.hg19)

# Load the data
path_pos = "C:\\Users\\khoah\\PhD_Documents\\STATS215\\final_project\\vdj_seqs_tsp7.fasta"
path_neg = "C:\\Users\\khoah\\PhD_Documents\\STATS215\\final_project\\non_vdj_seqs_tsp7.fasta"

seqs_pos = readDNAStringSet(path_pos)
seqs_neg = readDNAStringSet(path_neg)
gadem_vdj = GADEM(seqs_pos,verbose=1,genome=Hsapiens)
gadem_non_vdj = GADEM(seqs_neg,verbose=1,genome=Hsapiens)
length(seqs)

# function to scan PWM through a sequence and calculate probability
scan_pwm = function(pwm, seq) {
  n = nchar(seq)
  p = ncol(pwm)
  scores = numeric(n-p+1)
  for (i in 1:(n-p+1)) {
    subseq = substr(seq, i, i+p-1)
    subseq_onehot = onehot(subseq)
    scores[i] = sum(pwm[,1:p] * subseq_onehot)
  }
  return(scores)
}

# write function convert sequence to onehot representation order A,C,G,T
onehot = function(seq) {
  n = nchar(seq)
  onehot = matrix(0, n, 4)
  for (i in 1:n) {
    base = substr(seq, i, i)
    if (base == "A") {
      onehot[i,1] = 1
    } else if (base == "C") {
      onehot[i,2] = 1
    } else if (base == "G") {
      onehot[i,3] = 1
    } else if (base == "T") {
      onehot[i,4] = 1
    }
  }
  return(t(onehot))
}


#example
vdj_motifs = getPWM(gadem_vdj)
vdj_motifs[[1]]
seq = seqs[1]
scan_pwm(vdj_motifs[[1]], seq)

# scan all sequences with all motifs
vdj_scores = matrix(0, length(seqs), length(vdj_motifs))
for (i in 1:length(seqs)) {
  seq = seqs[i]
  for (j in 1:length(vdj_motifs)) {
    vdj_scores[i,j] = max(scan_pwm(vdj_motifs[[j]], seq))
  }
}
