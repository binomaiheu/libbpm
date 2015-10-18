
header-links: remove-links
	HEADERLIST="$(top_srcdir)/bpm*.h $(top_srcdir)/*/bpm*.h"; \
	for h in $$HEADERLIST; do \
	  BASENAME=`basename $$h`; \
	  test -r $$BASENAME || $(LN_S) $$h $$BASENAME; \
	done

remove-links: 
	rm -f bpm*.h


all: all-am header-links

clean: clean-am remove-links
distclean: distclean-am remove-links
	-rm -f Makefile
