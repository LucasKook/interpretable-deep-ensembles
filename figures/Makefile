
RFILES = $(wildcard *.R)

plots:
	for file in $(RFILES) ; do \
		R --vanilla < "$$file" ; \
	done

clean:
	rm *.pdf

all: plots clean
