R=R

all:bioqc-signedGenesets.html bioqc.html

clean:
	rm -f *.html
	rm -f *.md
	rm -rf *_files
	rm -rf *_cache

bioqc-signedGenesets.html:bioqc-signedGenesets.Rmd
	Rscript -e "rmarkdown::render('bioqc-signedGenesets.Rmd')"

bioqc.html:bioqc.Rmd
	Rscript -e "rmarkdown::render('bioqc.Rmd')"

