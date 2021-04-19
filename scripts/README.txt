###   Phylogenetic tree for Figure 2
Rscript Figure2.R

###	Boxplot for Figure 3a
Usage: Rscript  Make_Core_Pan.R InputFile OutputFile
Rscript Make_Core_Pan.R Figure3a.Data.txt Figure3a.svg

###	Enrichment test for NLR group per 1Mb window
Usage: perl EnrichedGroup.pl InputFolder/ Window Cutoff
perl EnrichedGroup.pl SupplementaryTable12.Data/ 1Mb 0.05

Note: 1) "Figure2.Data/" contains files of the phylogenetic tree (nwk format) and group information (csv format).
      2) "SupplementaryTable12.Data/" contains files of the number of NLRs assigned per 1 Mb window of the reference genome (CM334) for each accession.
      3) The output file names are "NLRsPerWindow.Total.1Mb.all.enriched.txt" and "NLRsPerWindow.Total.1Mb.enriched.txt".

###	Chi-square test for functional paralogs number variation (FPNV) of NLRs
Usage: Rscript Chi_square.FDR.R InputFile OutputFile Cutoff_of_standard_deviation
Rscript Chi_square.FDR.R SupplementaryTable13.Data.txt SupplementaryTable13.result.txt 2
