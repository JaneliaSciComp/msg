CFLAGS ?= -O2 -lm

.PHONY: all,msg,samtools,stampy,rpkgs
all: msg samtools stampy rpkgs

msg: countalleles hmmprobs

countalleles: countalleles.c
	$(CC) $(CFLAGS) -o $@ $<

hmmprobs: hmmprobs.c
	$(CC) $(CFLAGS) -o $@ $<

samtools: dependencies/samtools-0.1.9.tar.bz2
	tar -xjf $<
	cd samtools-0.1.9 && $(MAKE) && printf "Please set samtools_path in your msg.cfg to " && pwd | tr -d "\n" && printf "/samtools\n"

stampy: dependencies/stampy-1.0.32.tgz
	tar -xzf $<
	cd stampy-1.0.32 && $(MAKE) && printf "Stampy is installed to " && pwd

rpkgs: test_dependencies_R.sh
	sh $<

check: test_dependencies.sh
	sh $<
