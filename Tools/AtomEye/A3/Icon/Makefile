########################################################################
# Generate AXIcon data structure that contains a database of XPM icons #
########################################################################

PNG = $(wildcard Png/*.png)
XPM = $(addprefix Xpm/, $(addsuffix .xpm, $(basename $(notdir $(PNG)))))
SHELL = sh

.PHONY: all  i icon  c clean

all:	icon

i:	icon
icon:	$(XPM)
	@echo 'static char *icon[]={' > icon.c; \
	for xpm in Xpm/*.xpm; do \
	 sed -e 's/^[^\"].*//g' $$xpm >> icon.c; \
	 echo ',' >> icon.c; \
	 echo $$xpm added; \
	done; \
	echo '};' >> icon.c; \
	ls -la icon.c

Xpm/%.xpm:	Png/%.png
	convert $< $@

c:      clean
clean:
	@echo
	-rm -f $(PROG); \
	rm -f .\#* \#*\#; \
	rm -f a.out; \
	rm -f *~; \
	rm -f *.o *.png *.jpg *.eps; \
	rm -f *_TEST
	@echo
