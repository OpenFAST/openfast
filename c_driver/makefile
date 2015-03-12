PLATFORM = $(shell uname -s)

ifeq ($(PLATFORM),Darwin)
CC_TOOLS  = clang
CCFLAGS = -g -O1 -fsanitize=address -fno-omit-frame-pointer
else
CC_TOOLS  = gcc
#CCFLAGS = -g -O1 -fsanitize=address -fno-omit-frame-pointer
CCFLAGS = -g -std=c99 -O1 
endif

OBJ = test.o

all : $(OBJ)
	$(CC) $(CCFLAGS) -o c_driver $(OBJ) -L. ../src/libmap-1.00.01.so

clean:		
	rm -f *.o c_driver *.*~ *.*# makefile~ *.mod *.valgrind


memcheck:	
	valgrind --tool=memcheck --leak-check=yes --show-reachable=yes --track-origins=yes ./c_driver
	#valgrind --tool=exp-sgcheck ./c_driver
