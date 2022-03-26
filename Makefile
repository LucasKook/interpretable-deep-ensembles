
PKG = etram
VER = 0.0-1

build: $(PKG)
	R CMD build $<

install: $(PKG)_$(VER).tar.gz
	R CMD INSTALL $<

clean: $(PKG)_$(VER).tar.gz
	rm $<

all: build install clean
