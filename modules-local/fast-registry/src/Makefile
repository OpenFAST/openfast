ifeq ($(OS),Windows_NT)
  ifeq ($(OSTYPE),cygwin)
    RM=rm -f
  else
    RM=del
  endif
else
  RM = rm -f
endif
.SUFFIXES: .c .o

# i586-mingw32msvc-gcc
#CC_TOOLS = i586-mingw32msvc-gcc
CC_TOOLS = gcc
DEST_DIR = ..
CFLAGS = #-ansi
LDFLAGS =
DEBUG = -g
OBJ = registry.o \
	my_strtok.o \
	reg_parse.o \
	data.o \
	type.o \
	misc.o \
	sym.o \
	symtab_gen.o \
	gen_module_files.o \
	gen_c_types.o

# marco's: all : $(OBJ)
$(DEST_DIR)/registry.exe : $(OBJ)
	$(CC_TOOLS) -o $(DEST_DIR)/registry.exe $(DEBUG) $(LDFLAGS) $(OBJ)

.c.o :
	$(CC_TOOLS) $(CFLAGS) -c $(DEBUG) $<

clean clena:
	$(RM) $(OBJ) gen_comms.c standard.o

superclean : clean
	$(RM) $(DEST_DIR)/registry.exe  Registry_tmp.*

# regenerate this list with "makedepend -Y *.c"

# DO NOT DELETE THIS LINE -- make depend depends on it.

gen_module_files.o: protos.h registry.h data.h FAST_preamble.h

data.o: registry.h protos.h data.h
gen_allocs.o: protos.h registry.h data.h
gen_args.o: protos.h registry.h data.h
gen_scalar_derefs.o: protos.h registry.h data.h
gen_config.o: protos.h registry.h data.h
gen_defs.o: protos.h registry.h data.h
gen_mod_state_descr.o: protos.h registry.h data.h
gen_model_data_ord.o: protos.h registry.h data.h
gen_scalar_indices.o: protos.h registry.h data.h
gen_wrf_io.o: protos.h registry.h data.h
misc.o: protos.h registry.h data.h
my_strtok.o: registry.h protos.h data.h
reg_parse.o: registry.h protos.h data.h
registry.o: protos.h registry.h data.h Template_data.c Template_registry.c 
sym.o: sym.h
type.o: registry.h protos.h data.h
gen_interp.o: registry.h protos.h data.h
gen_streams.o: registry.h protos.h data.h
gen_c_types.o: Template_c_Types.c Template_c2f_helpers.c
