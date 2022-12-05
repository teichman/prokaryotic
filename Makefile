
# cmake is sufficiently annoying that we are wrapping it in Make

.PHONY: prokaryotic  # Always do cmake && make.
prokaryotic: # proto/addressbook.pb.h proto/addressbook.pb.cc
	mkdir -p build
	cd build && cmake ../ && make

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
	build/run

# proto/addressbook.pb.h proto/addressbook.pb.cc &: proto/addressbook.proto
# 	/opt/homebrew/Cellar/protobuf/21.9_1/bin/protoc -I=proto/ --cpp_out=proto/ proto/addressbook.proto

