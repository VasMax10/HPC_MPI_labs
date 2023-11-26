# size = 100

all: clean build gen run
std: build run

build: 
	g++ data_generation.cpp -o data_generation
	mpic++ parallel.cpp -o parallel
	g++ serial.cpp -o serial

gen:
	./data_generation $(size)

run:
	./serial $(size)
	mpirun --oversubscribe -np 2 ./parallel $(size) 
	mpirun --oversubscribe -np 4 ./parallel $(size) 

clean:
	rm -f matrix$(size).txt vector$(size).txt