SHELL := /bin/bash

download: GSE117715/GSE117715_Cynomolgus_monkey_aging_artery_count.txt.gz

samples: GSE117715/samples.csv

GSE117715/GSE117715_Cynomolgus_monkey_aging_artery_count.txt.gz: GSE117715
	wget -O $@ https://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117715/suppl/GSE117715_Cynomolgus_monkey_aging_artery_count.txt.gz

GSE117715/samples.csv: GSE117715/CellType_info.csv
	tail -n +2 $< | cut -d"," -f 2 | sort -u > $@

GSE117715:
	test -d $@ || mkdir $@
