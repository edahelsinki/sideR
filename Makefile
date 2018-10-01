FILES := Makefile R/run.R \
	R/sideR/global.R R/sideR/server.R R/sideR/ui.R R/sideR/sideR.R \
	R/runtime/runtime.R R/runtime/report.R \
	R/toyexample/toyexample.R \
	README.md LICENSE \
	R/sideR/data/toydata.rds \
	R/sideR/data/toy3.rds \
	R/sideR/data/segment.rds \
	R/sideR/data/bnc.rds \

dist/sideR.zip: $(FILES)
	rm -f dist/sideR.zip && zip -q dist/sideR.zip $(FILES)
