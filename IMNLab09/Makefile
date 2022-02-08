CXX = g++
EXEC = main.exe
CXXFLAGS = -std=c++17 -Ofast -lgsl -lgslcblas -lm
LIBS = -lm

SRC=$(wildcard *.cpp)
OBJS=$(SRC:.cpp=.o)
DEPS=$(SRC:.cpp=.d)

$(EXEC): $(OBJS) 
	$(CXX) $(CXXFLAGS) $+ -o $@


*.o: *.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean run rebuild

rebuild: 
	make clean
	make run

clean:
	rm -f *.o
	rm -f *~
	rm -f *.d
	rm -f *.dat
	rm -f *.png
	rm -f $(EXEC)

run: $(EXEC)
	./$(EXEC)
	gnuplot plot.sh

-include $(DEPS)
