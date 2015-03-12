CC  = g++
#NOTE: adrien @ marco : these are my release compile flags. You might want to change this
#CXXFLAGS = -O1 -fsanitize=address -fno-omit-frame-pointer -DNDEBUG -fmax-errors=5
CXXFLAGS = -O3 -DNDEBUG -fmax-errors=5
LDFLAGS = -L. ../src/libmap-1.00.01.so


OBJ = test.o

all : $(OBJ)
	$(CC) $(CXXFLAGS) -o cpp_driver $(OBJ) $(LDFLAGS)

clean:		
	rm -f *.o *.swp cpp_driver *.*~ *.*# makefile~ *.mod *.valgrind

memcheck:	
	valgrind --tool=memcheck --leak-check=yes --track-origins=yes ./cpp_driver
	#valgrind --tool=exp-sgcheck ./cpp_driver
