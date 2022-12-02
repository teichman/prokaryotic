prokaryotic: prokaryotic.cpp
	g++ -g -std=c++20 $^ \
	-I /opt/homebrew/Cellar/eigen/3.4.0_1/include/eigen3/ \
	-o $@

test: scratch_doctest
	./scratch_doctest

scratch_doctest: scratch_doctest.cpp
	g++ -g -std=c++20 $^ -o $@

clean:
	rm -f prokaryotic scratch_doctest
