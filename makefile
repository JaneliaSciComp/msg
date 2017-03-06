.PHONY: all,msg,samtools,rpkgs
all: msg samtools rpkgs

msg: countalleles hmmprobs

countalleles: countalleles.c
	$(CC) -O2 -lm -o $@ $<

hmmprobs: hmmprobs.c
	$(CC) -O2 -lm -o $@ $<

samtools: dependencies/samtools-0.1.9.tar.bz2
	tar -xjf $<
	cd samtools-0.1.9 && $(MAKE) && printf "Please set samtools_path in your msg.cfg to " && pwd

rpkgs: test_dependencies_R.sh
	sh $<

check: test_dependencies.sh
	sh $<