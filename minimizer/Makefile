CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++11 -I./include -Wno-sign-compare -Wno-deprecated-copy-with-user-provided-copy -Wno-unused-parameter
OBJ = ./bin/main.o ./bin/ref_minimizer.o ./bin/thr_minimizer.o ./bin/thr_minimizer_v2.o

default: ./bin/minimize_ref ./bin/minimize_thr_v1 ./bin/minimize_thr_v2
#	./$(TARGET)

./bin/minimize_ref: ./bin/main_ref.o ./bin/ref_minimizer.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_ref ./bin/main_ref.o ./bin/ref_minimizer.o

./bin/minimize_thr_v1: ./bin/main_thr_v1.o ./bin/thr_minimizer_v1.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_thr_v1 ./bin/main_thr_v1.o ./bin/thr_minimizer_v1.o

./bin/minimize_thr_v2: ./bin/main_thr_v2.o ./bin/thr_minimizer_v2.o
	$(CXX) $(CXXFLAGS) -o ./bin/minimize_thr_v2 ./bin/main_thr_v2.o ./bin/thr_minimizer_v2.o

./bin/main_ref.o: ./test/main_ref.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include

./bin/main_thr_v1.o: ./test/main_thr_v1.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include

./bin/main_thr_v2.o: ./test/main_thr_v2.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include

./bin/thr_minimizer_v1.o: ./c_code/thr_minimizer_v1.cpp ./c_code/functions.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include

./bin/thr_minimizer_v2.o: ./c_code/thr_minimizer_v2.cpp ./c_code/functions.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include

./bin/ref_minimizer.o: ./c_code/ref_minimizer.cpp ./c_code/functions.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@ -I./include


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
	rm -f $(OBJ) $(TARGET)

