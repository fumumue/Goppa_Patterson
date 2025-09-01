all:
	gcc -Wall -g -pg -O3 -mavx -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp Patterson.c

clang:
	clang -Wall -g -pg -O3 -mavx2 -mtune=native -march=native -ffast-math -funroll-loops  -fopenmp Patterson.c

clean:
	rm -f a.out

