BIOQC_DIR=/apps64/bi/ribios/BioQC/vignettes
BIOQC_EXAMPLE_DIR=/homebasel/biocomp/zhangj83/projects/2016-02-BioQCmanuscript/BioQC-example

all:example-kidney example-simulation bioqc bioqc-signedGenesets

example-kidney: ${BIOQC_EXAMPLE_DIR}/bioqc-kidney.md
	./fetch_content.sh $< 

example-simulation: ${BIOQC_EXAMPLE_DIR}/bioqc-simulation.md
	./fetch_content.sh $<

bioqc: ${BIOQC_DIR}/bioqc.md
	./fetch_content.sh $<

bioqc-signedGenesets: ${BIOQC_DIR}/bioqc-signedGenesets.md
	./fetch_content.sh $<

clean:
	rm -rfv ./pages/bioqc/*
