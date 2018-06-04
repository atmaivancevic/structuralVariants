#!/bin/bash

## Score SVs for potential pathogenicity
## Date: 4 June 2018 
##
## Example usage: ./scoreSVs.sh
## You will need the latest version of SVScore: https://github.com/lganel/SVScore

# remove chr perfix (for SVscore, which cannot have prefix)
cat germline.vcf | sed 's/chr//g' > germline_noPrefix.vcf

# use SVscore to add scores and gene exon/intron overlaps
svscore.pl -dv -e ~/Documents/Atma/references/genes/refGene.exons.bed -f ~/Documents/Atma/references/genes/refGene.introns.bed -c ~/Documents/Atma/references/snvs/whole_genome_SNVs.tsv.gz -i germline_noPrefix.vcf > germline_noPrefix_scored.vcf
