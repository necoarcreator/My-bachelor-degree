CC = g++
CFLAGS = -Wall 
All: execute

execute: main.o FT.o  
	$(CC) $(CFLAGS) main.o FT.o -o execute

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp -o main.o

FT.o: FT.cpp
	$(CC) $(CFLAGS) -c FT.cpp -o FT.o

clean: 
	@echo "Cleaning up..."
	rm execute
	rm -rf *.o execute 
	

rebuild: clean execute

run: execute
	./execute
