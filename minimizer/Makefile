CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++17 -I./include -I./src -I./test -I../tools/src/files -Wno-sign-compare -Wno-deprecated-copy-with-user-provided-copy -Wno-unused-parameter
LDFLAGS = -lz -pthread
OBJ = ./bin/main.o ./bin/ref_minimizer.o ./bin/thr_minimizer.o ./bin/thr_minimizer_v2.o ./bin/thr_minimizer_v3.o

# Create bin directory if it doesn't exist
$(shell mkdir -p ./bin)

default: ./bin/minimize_ref ./bin/minimize_thr_v1 ./bin/minimize_thr_v2 ./bin/minimizer_thr_v3
#	./$(TARGET)

./bin/minimize_ref: ./bin/main_ref.o ./bin/ref_minimizer.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_ref ./bin/main_ref.o ./bin/ref_minimizer.o $(LDFLAGS)

./bin/minimize_thr_v1: ./bin/main_thr_v1.o ./bin/thr_minimizer_v1.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_thr_v1 ./bin/main_thr_v1.o ./bin/thr_minimizer_v1.o $(LDFLAGS)

./bin/minimize_thr_v2: ./bin/main_thr_v2.o ./bin/thr_minimizer_v2.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_thr_v2 ./bin/main_thr_v2.o ./bin/thr_minimizer_v2.o $(LDFLAGS)

./bin/main_ref.o: ./test/main_ref.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/main_thr_v1.o: ./test/main_thr_v1.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/main_thr_v2.o: ./test/main_thr_v2.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/minimizer_thr_v3: ./bin/main_thr_v3.o ./bin/thr_minimizer_v3.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimizer_thr_v3 ./bin/main_thr_v3.o ./bin/thr_minimizer_v3.o $(LDFLAGS)

./bin/main_thr_v3.o: ./test/main_thr_v3.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/thr_minimizer_v1.o: ./src/thr_minimizer_v1.cpp ./src/functions.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/thr_minimizer_v2.o: ./src/thr_minimizer_v2.cpp ./src/functions.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/thr_minimizer_v3.o: ./src/thr_minimizer_v3.cpp ./src/functions.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

./bin/ref_minimizer.o: ./src/ref_minimizer.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


#
#
#

./bin/minimizer_hls_v2.o: ./hls/minimizer_v2.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include

./bin/minimize_hls_v2: ./bin/minimizer_hls_v2.o ./bin/main_hls_v2.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_thr_v2 ./bin/main_thr_v2.o ./bin/thr_minimizer_v2.o

#
#
#

clean:
	rm -rf ./bin/*.o ./bin/minimize_* ./bin/minimizer_*
