#!/bin/bash

# Function to extend peaks and get fasta
# Usage: extendPeaksAndGetFasta <peaks.narrowPeak.gz, chr, start, end, genome.fasta>
# Description: this function takes the narrowPeak as input, filters peaks in a region of interest, extends each peak by 100 bp each side and extract the sequence of each peak. The output files (bed and fasta) are written to the current directory.
# Main output: <peaks_extended.fasta>

extendPeaksAndGetFasta() {
  # Peaks in the region are selected and extended 100bp each direction
  gunzip -c $1 | awk -v chr=$2 -v start=$3 -v end=$4 -v OFS="\t" '$1==chr && $2 < end && $3 > start {print $1, $2 - 100, $3 + 100, $4}' > $(basename ${1/.narrowPeak.gz/_extended.bed})
  echo "Extended bed is in $PWD/$(basename ${1/.narrowPeak.gz/_extended.bed})"
  # The fasta is extracted
  bedtools getfasta -fi $5 -bed $(basename ${1/.narrowPeak.gz/_extended.bed}) > $(basename ${1/.narrowPeak.gz/_extended.fa})
  echo "Fasta is in $PWD/$(basename ${1/.narrowPeak.gz/_extended.fa})"
}

# Function to reformat the CTCF motif table copied from http://insulatordb.uthsc.edu (without header)
# Usage: reformat_insulatordb_table <*_extended_insdb_output.txt> <*_noEMBL.bed>
# Description: this function takes the insdb table as input, removes the EMBL motives (because very short), put . if motif score is negative and keeps the best motif for each peak.
# Output: <*_noEMBL.bed>

reformat_insulatordb_table() {
    grep -v EMBL "$1" | awk -v OFS="\t" '
{
  # $3 contains the coordinates of the peak
  split($3, a, ":|-")
  cur_chr = a[1]
  # cur_start contains the middle of the motif
  cur_start = a[2] + $4 + $5/2
  cur_motif = $1
  cur_score = $7
  cur_strand = $6
  # For each peak we only store the best motif
  if ($3 in scores_per_region){
    if (scores_per_region[$3] < cur_score){
      scores_per_region[$3] = cur_score
      infos_per_region[$3]["chr"] = cur_chr
      infos_per_region[$3]["start"] = cur_start
      infos_per_region[$3]["motif"] = cur_motif
      infos_per_region[$3]["score"] = cur_score
      infos_per_region[$3]["strand"] = cur_strand
    }
  } else {
    scores_per_region[$3] = cur_score
    infos_per_region[$3]["chr"] = cur_chr
    infos_per_region[$3]["start"] = cur_start
    infos_per_region[$3]["motif"] = cur_motif
    infos_per_region[$3]["score"] = cur_score
    infos_per_region[$3]["strand"] = cur_strand
  }
}
END{
  for (region in infos_per_region){
    # Multiple peaks may overlay, what is common is the motif center
    # If both peaks have the same motif we report it only once:
    name = infos_per_region[region]["chr"]"_"infos_per_region[region]["start"]
    if (!(name in written)){
      if (infos_per_region[region]["score"] < 0){
        # If the score is negative no orientation is given
        cur_strand = "."
      } else {
        cur_strand = infos_per_region[region]["strand"]
      }
      printf "%s\t%d\t%d\t%s\t%f\t%s\n",
        infos_per_region[region]["chr"], \
        infos_per_region[region]["start"], \
        infos_per_region[region]["start"] + 1, \
        region"_"infos_per_region[region]["motif"], \
        infos_per_region[region]["score"], cur_strand
      written[name] = 1
    }
  }
}' | sort -k1,1 -k2,3n > "$2"
}

# Function to merge replicates by keeping only the motifs that are present in all replicates
# Usage: keepMotivesRepIntersection <cat_sorted_allPeaks_allReps.bed> <nbReps> <cat_sorted_allReps.bed>
# Description: this function takes the concatenate file of all replicates and the number of replicates as input. It keeps only motives found in all replicates, and keeps the one with the best score.
# Output: <cat_sorted_allReps.bed>

keepMotivesRepIntersection() {
  awk -v nrep=$2 '
  {
    if(cur_chr == $1 && cur_start >= $2 - 1 && cur_start <= $2 + 1 && cur_strand == $6){
      timesFound=timesFound + 1
      if(cur_score < $5){
        bestLine = $0
        cur_score = $5
      }
    } else {
      if (NR > 1 && timesFound == nrep){
        print bestLine
      }
      bestLine = $0
      cur_chr = $1
      cur_start = $2
      cur_score = $5
      cur_strand = $6
      timesFound = 1
    }
  }' "$1" > "$3"
}

# Function to merge replicates by keeping the motifs that are present in at least one replicate
# Usage: keepMotivesRepUnion <cat_sorted_allPeaks_allReps.bed> <cat_sorted_atLeastOneRep.bed>
# Description: this function takes the concatenate file of all replicates. For each motif, it keeps the one with the best score.
# Output: <cat_sorted_atLeastOneRep.bed>

keepMotivesRepUnion() {
  awk '
  {
    if(cur_chr == $1 && cur_start >= $2 - 1 && cur_start <= $2 + 1 && cur_strand == $6){
      if(cur_score < $5){
        bestLine = $0
        cur_score = $5
      }
    } else {
      if (NR > 1){
        print bestLine
      }
      bestLine = $0
      cur_chr = $1
      cur_start = $2
      cur_score = $5
      cur_strand = $6
    }
  }' "$1" > "$2"
}

# Function to create a bed with rgb field corresponding to motif orientation of CTCF
# Usage: makeRgbField <*extended_noEMBL.bed> [inverted]
# Description: this function takes the bed6 as input, and return a bed9 as output with a different color for positive and negative motives. Positive motives are red (236,28,36), Negative motives are blue (46,49,145) except if 'inverted' is set, then it is the contrary. The output is written to the current directory.
# Output: <*extended_noEMBL_colored.bed>

makeRgbField() {
  awk -F "\t" -v OFS="\t" -v xaxis="$2" -v colorRed="236,28,36" -v colorBlue="46,49,145" -v colorOther="0,0,0" '
  {
    if ($6 == "+"){
      color = colorRed
      if (xaxis == "inverted"){
        color = colorBlue
      }
    } else if ($6 == "-" ){
      color = colorBlue
      if (xaxis == "inverted"){
        color = colorRed
      }
    } else {
      color = colorOther
    }
    print $1, $2, $3, $4, $5, $6, $2, $2, color
  }' "$1" > $(basename "${1/.bed/_colored.bed}")
}
