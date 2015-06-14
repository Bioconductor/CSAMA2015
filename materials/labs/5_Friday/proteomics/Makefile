proteomics.pdf: proteomics.Rmd
	Rscript -e 'rmarkdown::render("proteomics.Rmd", output_format="all")'
	Rscript -e 'knitr::purl("proteomics.Rmd")'

clean:
	rm -f erwinia_carotovora.fasta
	rm -f erwinia_carotovora.revCat.*
	rm -f TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.*
	rm -f F063721.dat-mztab.txt
	rm -rf proteomics_files proteomics_cache .Rcache

