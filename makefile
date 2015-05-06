all:
	g++ -lm nxn.cpp -o output

run:
	./output > output.txt
