#
default: install
#
OBJFILES= m_dom_error.$(OBJEXT) m_dom_parse.$(OBJEXT) m_dom_utils.$(OBJEXT) \
	  m_dom_extras.$(OBJEXT) m_dom_dom.$(OBJEXT) FoX_dom.$(OBJEXT)
#
LIBRARY=libFoX_dom.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
include ../arch.make
#
#
m_dom_dom.F90: $(M4FILES)
#
CP=cp -p
install: $(OBJFILES)
	@echo "  ==> Creating $(LIBRARY) with $(OBJFILES)"
	$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR)
	for i in *.$(MOD_EXT) ; do  \
           $(CP) $$i $(MOD_DIR); done
#
check:
	
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
# DO NOT DELETE THIS LINE - used by make depend
FoX_dom.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT) m_dom_extras.$(OBJEXT) m_dom_parse.$(OBJEXT) m_dom_utils.$(OBJEXT)
m_dom_dom.$(OBJEXT): m_dom_error.$(OBJEXT)
m_dom_extras.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT)
m_dom_parse.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT)
m_dom_utils.$(OBJEXT): m_dom_dom.$(OBJEXT) m_dom_error.$(OBJEXT)
