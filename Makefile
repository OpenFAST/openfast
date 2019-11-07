MAKE=make --no-print-directory
MAKE_ARGS=VERBOSE=1
CMAKE_ARGS=-DCMAKE_BUILD_TYPE=DEBUG -DDOUBLE_PRECISION:BOOL=OFF -DGENERATE_TYPES:BOOL=ON  
suffix=-single-db
ifeq ($(OS),Windows_NT)
	OS=Windows
	suffix=-nmake$(suffix)
	CMAKE_ARGS=$(CMAKE_ARGS) -G"NMake Makefiles"
	MAKE=nmake
	RMDIR=rmdir /S /Q
	LIBEXT='.dll'
	EXEEXT='.exe'
else
    OS=$(shell uname -s)
    ifeq ($(OS),Linux)
		# ISSUE WITH CMAKE PICKING UP F95 TODO
		FC=/usr/bin/gfortran
# 		FC=/usr/bin/gfortran-7
		LIBEXT='.so'
		LIBPRE='lib'
    else ifeq ($(OS),Darwin)
		LIBEXT='.dylib'
    endif
	RMDIR=rm -rf
endif
BUILD_DIR=build$(suffix)
TEST_DIR=../tcf-test


all: compile copy test


$(BUILD_DIR):
	@echo "------------------------------------------------------------"
	mkdir $(BUILD_DIR)

compile: $(BUILD_DIR)
	@echo "------------------------------------------------------------"
	cd $(BUILD_DIR) && cmake $(CMAKE_ARGS) .. && $(MAKE) $(MAKE_ARGS)

clean:
	cd $(BUILD_DIR) && $(MAKE) clean

copy:
	cp $(BUILD_DIR)/glue-codes/openfast/openfast$(EXEEXT) $(TEST_DIR)

test:
	cd $(TEST_DIR) && $(MAKE)
	
