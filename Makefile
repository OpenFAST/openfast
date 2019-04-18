MAKE=make --no-print-directory
ifeq ($(OS),Windows_NT)
	OS=Windows
	suffix=-nmake
	CMAKE_ARGS=-G"NMake Makefiles"
	MAKE=nmake
	RMDIR=rmdir /S /Q
	LIBEXT='.dll'
	EXEEXT='.exe'
else
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
endif
BUILD_DIR=build$(suffix)
TEST_DIR=../noise-test/



all: compile copy test


$(BUILD_DIR):
	@echo "------------------------------------------------------------"
	mkdir $(BUILD_DIR)

compile: $(BUILD_DIR)
	@echo "------------------------------------------------------------"
	cd $(BUILD_DIR) && cmake $(CMAKE_ARGS) .. && $(MAKE)

clean:
	cd $(BUILD_DIR) && $(MAKE) clean

copy:
	cp $(BUILD_DIR)/glue-codes/openfast/openfast$(EXEEXT) $(TEST_DIR)

test:
	cd $(TEST_DIR) && $(MAKE)
	

debug:
	build/glue-codes/openfast/openfast SC_FAST_A001/OC4Jacket.fst > OUT_NEW
	build/glue-codes/openfast/openfast SC_FAST_A001/OC4Jacket_newref.fst > OUT

