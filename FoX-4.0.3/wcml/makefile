#
default: install
#
OBJFILES=FoX_wcml.$(OBJEXT) m_wcml_coma.$(OBJEXT) m_wcml_core.$(OBJEXT) m_wcml_stml.$(OBJEXT) \
         m_wcml_parameter.$(OBJEXT) m_wcml_property.$(OBJEXT) m_wcml_metadata.$(OBJEXT) \
         m_wcml_lattice.$(OBJEXT) m_wcml_geometry.$(OBJEXT) m_wcml_molecule.$(OBJEXT) \
         m_wcml_lists.$(OBJEXT)
#
LIBRARY=libFoX_wcml.$(LIBEXT)
INCFLAGS=$(MOD_PREFIX)../objs/finclude
#
#
include ../arch.make
#
CP=cp -p
install: $(OBJFILES)
	@echo "  ==> Creating $(LIBRARY) with $(OBJFILES)"
	@$(AR) $(ARFLAGS_EXTRA) cru $(LIBRARY) $(OBJFILES)
	@$(RANLIB) $(LIBRARY)
	$(CP) $(LIBRARY) $(LIB_DIR) ;
	for i in *.$(MOD_EXT) ; do  \
           $(CP) $$i $(MOD_DIR); done
#
check:
	(cd test; make) | tee -a ../check.out
#
clean:
	rm -f *.$(OBJEXT)  *.$(MOD_EXT) $(LIBRARY)
#
# DO NOT DELETE THIS LINE - used by make depend
FoX_wcml.$(OBJEXT): m_wcml_coma.$(OBJEXT) m_wcml_core.$(OBJEXT) m_wcml_geometry.$(OBJEXT) m_wcml_lattice.$(OBJEXT)
FoX_wcml.$(OBJEXT): m_wcml_lists.$(OBJEXT) m_wcml_metadata.$(OBJEXT) m_wcml_molecule.$(OBJEXT)
FoX_wcml.$(OBJEXT): m_wcml_parameter.$(OBJEXT) m_wcml_property.$(OBJEXT) m_wcml_stml.$(OBJEXT)
m_wcml_coma.$(OBJEXT): m_wcml_stml.$(OBJEXT)
m_wcml_core.$(OBJEXT): m_wcml_metadata.$(OBJEXT)
m_wcml_parameter.$(OBJEXT): m_wcml_stml.$(OBJEXT)
m_wcml_property.$(OBJEXT): m_wcml_stml.$(OBJEXT)
