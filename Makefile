CXX = g++
CXXFLAGS = -Wall -Wextra -std=c++11
LIBS = -larmadillo
LDFLAGS =
TARGET = bin/exec
OBJS = obj/main.o obj/Poly.o obj/basis.o
CXXTEST = cxxtestgen
CXXTESTFLAGS = --error-printer
TARGET_TEST = bin/test
TESTS = tests/MandatoryTests.h
OBJS_TEST = tests/tests.o obj/Poly.o obj/basis.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

obj/%.o: src/%.cpp headers/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

obj/main.o: src/main.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


#Tests

test: $(TARGET_TEST)
	rm -f bin/test.cpp

$(TARGET_TEST): $(OBJS_TEST)
	$(CXX) $^ -o $@ $(LIBS) 

tests/test.o: bin/test.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

tests/%.cpp: $(TESTS)
	$(CXXTEST) $(CXXTESTFLAGS) -o $@ $<


.PHONY: clean all
clean:
	rm -f ./obj/*.o
	rm -f ./filesResult/*.csv
	rm -f ./tests/*.o
	rm -f ./bin/*
	rm -f sphere.df3
	rm -f visu.png
	rm -f result.png
	rm -f result.csv