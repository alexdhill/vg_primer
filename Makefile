GPP = g++ -g -Wall -std=c++14

all: filter_primers

clean:
	rm *.o filter_primers

vgfilters.o: vgfilters.cpp vgfilters.hpp
	$(GPP) -c vgfilters.cpp -o vgfilters.o

read_primer3.o: read_primer3.cpp read_primer3.hpp
	$(GPP) -c read_primer3.cpp -o read_primer3.o

main.o: main.cpp read_primer3.o
	$(GPP) -c main.cpp -o main.o

filter_primers: main.o read_primer3.o vgfilters.o
	$(GPP) -o filter_primers main.o read_primer3.o vgfilters.o -lsdsl -lvgio -lbdsg -lhandlegraph -lvcflib -ltabixpp -lgssw -lssw -lsublinearLS -lpthread -lncurses -lgcsa2 -lgbwtgraph -lgbwt -lkff -ldivsufsort -ldivsufsort64 -lvcfh -lraptor2 -lpinchesandcacti -l3edgeconnected -lsonlib -lfml -lstructures -lxg -lzstd -lgtest -lprotobuf -ljansson -lomp; \
	rm *.o