SHELL := /bin/bash

# Get Samples
SYAML = $(shell find ../data/SampleMetadata -name "*.yaml")
SID = $(basename $(notdir $(SYAML)))
SIDX1 = $(SID:%=Samples/CA_%.yaml)
SIDX2 = $(SID:%=Samples/AA_%.yaml)

# Get Contrasts
CONT = $(shell cat ../settings/contrasts.txt)
CIDX = $(CONT:%=Contrasts/%.yaml)

# Get Contrasts
GR = $(shell cat ../settings/groups.txt)
GIDX = $(GR:%=Groups/%.yaml)

.PHONY: samples groups contrasts

# Index samples

samples: samples1 samples2

samples1: $(SIDX1)

Samples/CA_%.yaml: ../data/SampleMetadata/%.yaml ../settings/QC.yaml ../settings/seurat.yaml Samples
	@echo Creating $@
	@grep -v "^\.\.\." $< > $@ 
	@grep -Ev "^(\.\.\.|---)" $(word 2, $^) >> $@ 
	@grep -Ev "^(\.\.\.|---)" $(word 3, $^) >> $@ 
	@echo "..." >> $@

samples2: $(SIDX2)

Samples/AA_%.yaml: ../data/SampleMetadata/%.yaml ../settings/QC.yaml ../settings/seurat.yaml Samples
	@echo Creating $@
	@grep -v "^\.\.\." $< > $@ 
	@grep -Ev "^(\.\.\.|---)" $(word 2, $^) >> $@ 
	@grep -Ev "^(\.\.\.|---)" $(word 3, $^) >> $@ 
	@echo "..." >> $@

Samples:
	test -d $@ || mkdir $@

# Index contrasts

contrasts: $(CIDX)

comma:=$(shell echo ", ")

Contrasts/%.yaml: ../settings/seurat.yaml Contrasts
	@echo Creating $@
	@echo "---" > $@
	@echo 'groups: [$(subst -,$(comma),$*)]' >> $@
	@echo 'annotation: Annotation.Level.1' >> $@
	@grep -Ev "^(\.\.\.|---)" $< >> $@ 
	@echo "..." >> $@

Contrasts:
	test -d $@ || mkdir $@

# Index groups 

groups: $(GIDX)

Groups/%.yaml: $(YAML) Groups
	@echo Creating $@
	@echo "---" > $@
	@echo 'samples: [$(shell grep -r -B 2 "biotype.*$*" ../data/SampleMetadata/ | grep "id" | sed -e "s/.*: *//g" | awk '{if(NR==1) c=$$0; else c=c", "$$0}END{print c;}')]' >> $@
	@echo 'annotation: Annotation.Level.1' >> $@
	@echo "..." >> $@

Groups:
	test -d $@ || mkdir $@
