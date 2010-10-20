default:

install-structures:
	@if [ "x${QUIP_STRUCTS_DIR}" == "x" ]; then \
	  if [ "x${QUIP_DIR}" == "x" ]; then \
	    if [ -d ${HOME}/share/quip_structures ]; then \
	      echo "QUIP_STRUCTURS_DIR, QUIP_DIR not defined, copying structures into HOME/share/quip_structures"; \
	      cp *xyz ${HOME}/share/quip_structures; \
	    else \
	      echo "No QUIP_STRUCTS_DIR or QUIP_DIR defined, and HOME/share/quip_structures is not a directory"; \
	    fi; \
	  else \
	    if [ -d ${QUIP_DIR}/structures ]; then \
	      echo "QUIP_STRUCTS_DIR not defined, copying structures into QUIP_DIR/structures"; \
	      cp *xyz ${QUIP_DIR}/structures; \
	    else \
	      echo "QUIP_DIR defined as ${QUIP_DIR}, but ${QUIP_DIR}/structures is not a directory"; \
	    fi; \
	  fi; \
	else \
	  if [ -d ${QUIP_STRUCTS_DIR} ]; then \
	    echo "copying structures into QUIP_STRUCTS_DIR"; \
	    mkdir -p ${QUIP_STRUCTS_DIR}/structures; \
	    cp *xyz ${QUIP_STRUCTS_DIR}/structures; \
	  else \
	    echo "QUIP_STRUCTS_DIR defined as ${QUIP_STRUCTS_DIR}, but it is not a directory"; \
	  fi; \
	fi
