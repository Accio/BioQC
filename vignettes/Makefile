R=R

all:bioqc-signedGenesets.html bioqc.html

clean:
	rm -fv *.html
	rm -fv *.md
	rm -rfv *_files
	rm -rfv *_cache

bioqc-signedGenesets.html:bioqc-signedGenesets.Rmd
	Rscript -e "rmarkdown::render('bioqc-signedGenesets.Rmd')"

bioqc.html:bioqc.Rmd
	Rscript -e "rmarkdown::render('bioqc.Rmd')"

