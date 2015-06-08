CFLAGS = -O2
VERSION = 1.1
BINDIR = /usr/bin
MANDIR = /usr/man/man1
PAGER = less

DIST= asa.c asa.1 Makefile COPYING asa.dat README

all: asa

install: all
	install -s -m 755 asa $(BINDIR)
	install -m 644 asa.1 $(MANDIR)

demo: all
	./asa asa.dat | $(PAGER)

asa.o: asa.c
	$(CC) $(CFLAGS) -c asa.c

asa: asa.o
	$(CC) $(CFLAGS) -o asa asa.o

checkin: $(DIST)
	ci -l $(DIST)

dist: $(DIST)
	(cd ..; tar cf - `for a in $(DIST); do echo asa-$(VERSION)/$$a; done`) |\
	gzip -9 > ../asa-$(VERSION).tar.gz
