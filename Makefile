
msSrc ?= ..


all:  ; 
	$(MAKE) -f  ${msSrc}/script/Makefile.in  recurseMake USE_msRecurse=1
	cp AllRunContAngle  ${msSrc}/../bin/foamx4m/

clean:; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseClean USE_msRecurse=1

tsts=test.sh
USE_msTEST=1
include  ${msSrc}/script/Makefile.in
