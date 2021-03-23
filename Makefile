
.PHONY: settings index

all: settings index

# data:
# 	$(MAKE) -C data download
# 	$(MAKE) -C data organize 
# 	$(MAKE) -C data clean 

settings:
	$(MAKE) -C settings zhang
	$(MAKE) -C settings groups.txt
	$(MAKE) -C settings contrasts.txt

index:
	$(MAKE) -C index samples 
	$(MAKE) -C index groups
	$(MAKE) -C index contrasts
