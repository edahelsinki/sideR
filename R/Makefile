FILES := Makefile run.R \
	sideR/global.R sideR/server.R sideR/ui.R sideR/sideR.R \
	runtime/runtime.R runtime/report.R \
	toyexample/toyexample.R \
	README.txt \
	sideR/data/toydata.rds \
	sideR/data/toy3.rds \
	sideR/data/segment.rds \
	sideR/data/bnc.rds \
	# sideR/data/data_122_tm_alig.rds \
	# sideR/data/data_211_dia_alig.rds \
	# sideR/data/celldata.rds \

sideR.zip: $(FILES)
	rm -f sideR.zip && zip -q sideR.zip $(FILES)
