all:
	R CMD build S4v2 --no-manual
	R CMD check S4v2 --no-manual
install:
	R CMD INSTALL S4v2
clean:
	rm -f S4v2/src/S4v2.o S4v2/src/S4v2.so S4v2_*.tar.gz
	rm -r S4v2.Rcheck