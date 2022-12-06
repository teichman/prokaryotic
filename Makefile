
# cmake is sufficiently annoying that we are wrapping it in Make

.PHONY: prokaryotic  # Always do cmake && make.
prokaryotic:
	mkdir -p build
	cd build && cmake ../ && make -j 8

clean:
	rm -rf build CMakeFiles
	rm -f proto/*.pb.h proto/*.pb.cc CMakeCache.txt

test: prokaryotic
	@echo
	@echo "============================================================"
	@echo "= Running tests"
	@echo "============================================================"
	build/test -s

run: prokaryotic
	@echo
	@echo "============================================================"
	@echo "= Running simulation"
	@echo "============================================================"
	build/run config.yaml

debug_snippet: prokaryotic
	lldb -o run -- build/debug_snippet

