ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ARCH = debug
NDIM := 3
NDIR := 3

export ARCH NDIM NDIR

SRC_DIR := $(AMRVAC_DIR)/src
LIB_DIR := $(AMRVAC_DIR)/lib/$(NDIM)d$(NDIR)_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make
LIB_AMRVAC := $(LIB_DIR)/libamrvac.a

INC_DIRS := $(LIB_DIR)
LIB_DIRS := $(LIB_DIR)
LIBS := amrvac

.PHONY: all clean allclean force hdr

all: | hdr amrvac

hdr:
ifeq ("$(wildcard amrvac.h)","")
	touch amrvac.h
endif
	@mkdir -p $(LIB_DIR)
	rsync -c amrvac.h $(LIB_DIR)/amrvac.h

include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

-include local.make

vpath %.t $(SRC_DIR)

.PRECIOUS: mod_usr.f
.INTERMEDIATE: amrvac.o mod_usr.o mod_usr.mod

$(LIB_AMRVAC): force
	@mkdir -p $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE)

clean:
	@echo 'Cleaning local objects ("make allclean" cleans libamrvac)'
	$(RM) amrvac amrvac.o mod_usr.o mod_usr.mod

allclean:
	$(RM) amrvac amrvac.o mod_usr.o mod_usr.f mod_usr.mod
	@echo 'Cleaning libamrvac'
	@mkdir -p $(LIB_DIR)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE) clean

amrvac: mod_usr.o amrvac.o

amrvac.o mod_usr.o: $(LIB_AMRVAC)
amrvac.o: mod_usr.o
mod_usr.o: mod_usr.f amrvac.h
