default:

install:
	@if [[ -z "${QUIP_STRUCTS_DIR}" ]]; then \
	  if [[ -d ${HOME}/share/quip_structures ]]; then \
	    echo "QUIP_STRUCTS_DIR not defined, copying structures into ${HOME}/share/quip_structures"; \
	    cp *xyz ${HOME}/share/quip_structures; \
	  else \
	    echo "QUIP_STRUCTS_DIR not defined and ${HOME}/share/quip_structures does not exist or is not a directory"; \
	  fi; \
	else \
	  if [[ -e "${QUIP_STRUCTS_DIR}" && ! -d "${QUIP_STRUCTS_DIR}" ]]; then \
	    echo "QUIP_STRUCTS_DIR defined as ${QUIP_STRUCTS_DIR} exists, but it is not a directory"; \
	  else \
	    echo "copying structures into QUIP_STRUCTS_DIR=${QUIP_STRUCTS_DIR}"; \
	    mkdir -p ${QUIP_STRUCTS_DIR}; \
	    cp *xyz ${QUIP_STRUCTS_DIR}; \
	  fi; \
	fi
