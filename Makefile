BIOQC_DIR=/homebasel/biocomp/sturmg/projects/BioQC/vignettes
BIOQC_EXAMPLE_DIR=/homebasel/biocomp/sturmg/projects/BioQC-example

all: examples vignettes

examples: example-kidney example-simulation example-wmw-performance
vignettes: bioqc bioqc-signedGenesets bioqc-efficiency

example-kidney: ${BIOQC_EXAMPLE_DIR}/bioqc-kidney.md
	./fetch_content.sh $< 

example-simulation: ${BIOQC_EXAMPLE_DIR}/bioqc-simulation.md
	./fetch_content.sh $<

example-wmw-performance: ${BIOQC_EXAMPLE_DIR}/bioqc-wmw-test-performance.md
	./fetch_content.sh $<

bioqc: ${BIOQC_DIR}/bioqc.md
	./fetch_content.sh $<

bioqc-signedGenesets: ${BIOQC_DIR}/bioqc-signedGenesets.md
	./fetch_content.sh $<

bioqc-efficiency: ${BIOQC_DIR}/bioqc-efficiency.md
	./fetch_content.sh $<

clean:
	rm -rfv ./pages/bioqc/*

serve:
	bundle exec jekyll serve
