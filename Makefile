
# cmake is sufficiently annoying that we are wrapping it in Make

.PHONY: rel  # Always do cmake && make.
rel:
	mkdir -p build/$@
	cd build/$@ && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../../ && make -j 8

.PHONY: debug  # Always do cmake && make.
debug:
	mkdir -p build/$@
	cd build/$@ && cmake -DCMAKE_BUILD_TYPE=Debug ../../ && make -j 8

clean:
	rm -rf build CMakeFiles
	rm -f proto/*.pb.h proto/*.pb.cc CMakeCache.txt

test: rel
	@echo
	@echo "============================================================"
	@echo "= Running tests in release mode"
	@echo "============================================================"
	build/rel/test -s

test-debug: debug
	@echo
	@echo "============================================================"
	@echo "= Running tests in debug mode"
	@echo "============================================================"
	build/debug/test -s

