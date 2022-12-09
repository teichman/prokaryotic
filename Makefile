
# cmake is sufficiently annoying that we are wrapping it in Make

.PHONY: build/debug  # Always do cmake && make.
build/debug:
	mkdir -p $@
	cd $@ && cmake -DCMAKE_BUILD_TYPE=Debug ../../ && make -j 8

.PHONY: build/rel  # Always do cmake && make.
build/rel:
	mkdir -p $@
	cd $@ && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ../../ && make -j 8

clean:
	rm -rf build CMakeFiles
	rm -f proto/*.pb.h proto/*.pb.cc CMakeCache.txt

test-rel: build/rel
	@echo
	@echo "============================================================"
	@echo "= Running tests in release mode"
	@echo "============================================================"
	build/rel/test -s

test-debug: build/debug
	@echo
	@echo "============================================================"
	@echo "= Running tests in debug mode"
	@echo "============================================================"
	build/debug/test -s

