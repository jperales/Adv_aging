SHELL := /bin/bash

# Get Samples
SYAML = $(shell find ../data/SampleMetadata -name "*.yaml")

clean:
	rm contrasts.txt groups.txt

contrasts.txt:
	@echo Creating '$@'
	@echo old-young > $@		

groups.txt: ../data/SampleMetadata $(SYAML)
	@echo Creating '$@'
	@grep -r "biotype" $< | sed -e "s/.*: *//g" | sort -u > $@ 
