default:

install-dtds:
	@if [[ -z "${QUIP_DTD_DIR}" ]]; then \
	  if [[ -d ${HOME}/share/quip_dtds ]]; then \
	    echo "QUIP_DTD_DIR not defined, copying dtds into HOME/share/quip_dtds"; \
	    cp *dtd ${HOME}/share/quip_dtds; \
	  else \
	    echo "QUIP_DTD_DIR not defined and HOME/share/quip_dtds does not exist or is not a directory"; \
	  fi; \
	else \
	  if [[ -e "${QUIP_DTD_DIR}" && ! -d "${QUIP_DTD_DIR}" ]]; then \
	    echo "QUIP_DTD_DIR defined as ${QUIP_DTD_DIR} exists, but it is not a directory"; \
	  else \
	    echo "copying dtds into QUIP_DTD_DIR=${QUIP_DTD_DIR}"; \
	    mkdir -p ${QUIP_DTD_DIR}; \
	    cp *dtd ${QUIP_DTD_DIR}; \
	  fi; \
	fi
