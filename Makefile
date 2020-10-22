CC = c++
CFLAGS = -mavx

all: main.cpp function.cpp utility.cpp
	$(CC) $(CFLAGS) -o run $^

# all: main.o function.o utility.o
# 	$(CC) $(CFLAGS) -o run $^

# main.o: main.cpp function.hpp
# function.o: function.cpp function.hpp
# 	$(CC) -c $^ -mavx
# utility.o: utility.cpp function.hpp
origin.o : origin.cpp function.hpp


.PHONY: clean
clean:
	# rm -f *.o 
	rm -f run