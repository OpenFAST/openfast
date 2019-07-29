MAKE=make --no-print-directory
ifeq ($(OS),Windows_NT)
	OS=Windows
	suffix=-nmake
	CMAKE_ARGS=-G"NMake Makefiles"
	MAKE=nmake
	RMDIR=rmdir /S /Q
	LIBEXT=.dll
	EXEEXT=.exe
	COPY=copy /y 
	SLASH=\\
else
    COPY=cp
    OS=$(shell uname -s)
	suffix=
    ifeq ($(OS),Linux)
		# ISSUE WITH CMAKE PICKING UP F95 TODO
		FC=/usr/bin/gfortran
		LIBEXT='.so'
		LIBPRE='lib'
    else ifeq ($(OS),Darwin)
		LIBEXT='.dylib'
    endif
	RMDIR=rm -rf
	SLASH="/"
endif
BUILD_DIR=build$(suffix)
TEST_DIR=..$(SLASH)openfast-noise-test$(SLASH)
CMAKE=cmake



all: compile copy test


$(BUILD_DIR):
	@echo "------------------------------------------------------------"
	mkdir $(BUILD_DIR)

compile: $(BUILD_DIR)
	@echo "------------------------------------------------------------"
	cd $(BUILD_DIR) && $(CMAKE) $(CMAKE_ARGS) .. && $(MAKE)

clean:
	cd $(BUILD_DIR) && $(MAKE) clean

copy:
	$(COPY) $(BUILD_DIR)$(SLASH)glue-codes$(SLASH)openfast$(SLASH)openfast$(EXEEXT) $(TEST_DIR)

test:
	cd $(TEST_DIR) && make
	

debug:
	build/glue-codes/openfast/openfast SC_FAST_A001/OC4Jacket.fst > OUT_NEW
	build/glue-codes/openfast/openfast SC_FAST_A001/OC4Jacket_newref.fst > OUT

