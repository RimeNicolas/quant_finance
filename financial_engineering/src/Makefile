HEADERS = motion.h random_generator.h print_func.h
OBJECTS = main.o 
SPATH = ./src/
EPATH = ./build/
OPATH = ./build/object/

default: main

%.o: ./%.cpp $(HEADERS)
	g++ -c $< -o $(OPATH)$@

main: $(OBJECTS)
	g++ $(OPATH)*.o -o $(EPATH)$@

clean:
	-rm -f $(OPATH)*.o
	-rm -f $(EPATH)main
