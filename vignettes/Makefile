R=R

all:bioqc-signedGenesets.html bioqc.html bioqc-efficiency.html

clean:
	rm -fv *.html
	rm -fv *.md
	rm -rfv *_files
	rm -rfv *_cache

bioqc-signedGenesets.html:bioqc-signedGenesets.Rmd
	Rscript -e "rmarkdown::render('bioqc-signedGenesets.Rmd', output_format='all')"

bioqc.html:bioqc.Rmd
	Rscript -e "rmarkdown::render('bioqc.Rmd', output_format='all')"

bioqc-efficiency.html:bioqc-efficiency.Rmd
	Rscript -e "rmarkdown::render('bioqc-efficiency.Rmd', output_format='all')"
