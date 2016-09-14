all:example-kidney example-simulation bioqc bioqc-signedGenesets

example-kidney: /homebasel/biocomp/sturmg/projects/BioQC-example/bioqc-kidney.md
	./fetch_content.sh $< 

example-simulation: /homebasel/biocomp/sturmg/projects/BioQC-example/bioqc-simulation.md
	./fetch_content.sh $<

bioqc: /homebasel/biocomp/sturmg/projects/BioQC/bioqc.md
	./fetch_content.sh $<

bioqc-signedGenesets: /homebasel/biocomp/sturmg/projects/BioQC/bioqc-signedGenesets.md
	./fetch_content.sh $<

clean:
	rm -rfv ./pages/bioqc/*
