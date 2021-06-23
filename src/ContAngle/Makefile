
msSrc ?= ..

all:    wmake


wmake:
	${msSrc}/script/AllMake
	cp AllRunContAngle  ${msSrc}/../bin/foamx3m/


clean:
	${msSrc}/script/AllClean


tsts=test.sh
USE_msTEST=1
include  ${msSrc}/script/Makefile.in
