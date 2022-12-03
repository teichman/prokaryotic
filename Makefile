
# cmake is sufficiently annoying that we are wrapping it in Make

.PHONY: prokaryotic  # Always do cmake && make.
prokaryotic:
	mkdir -p build
	cd build && cmake ../ && make

clean:
	rm -rf build

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
	build/run

