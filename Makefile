
.PHONY: data settings index

all: data settings index

data:
	$(MAKE) -C data download
	$(MAKE) -C data samples 

settings:
	$(MAKE) -C settings groups.txt
	$(MAKE) -C settings contrasts.txt

index:
	$(MAKE) -C index samples 
	$(MAKE) -C index groups
	$(MAKE) -C index contrasts
