prokaryotic: prokaryotic.cpp
	g++ -g -std=c++20 $^ \
	-I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ \
	-o $@



clean:
	rm prokaryotic
