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
	mpirun -n 4 ./parallel $(size) 

clean:
	rm -f AMatrix$(size).txt BMatrix$(size).txt