headers = ./header/*.h
files = ./*cpp ./test/*.cpp ./definition/*.cxx
objects = ./lib/Nodes.a ./lib/NodesHelper.a
Ldir = -L/usr/lib/x86_64-linux-gnu
libs = -lboost_regex
def = ./definition/
G++ = g++ -std=c++1y -Og -Wall -pthread -Wl,--no-as-needed $(Ldir)
Gdbg = g++ -std=c++1y -Og -Wall -Wshadow -ggdb -pthread -Wl,--no-as-needed $(Ldir)
libdir = lib

bindir = bin
release = $(bindir)/main.out

all: $(release)

gif: $(release)
	./$(release)
	./tools/make_gif.bash

$(release): $(files) $(objects) main.cpp $(bindir)
	$(G++) $(libs) $(objects) -o $(release) main.cpp

./lib/Nodes.a: $(headers) $(def)/Nodes.cxx $(libdir)
	$(G++) -o ./$(libdir)/Nodes.a -c $(def)/Nodes.cxx

./lib/NodesHelper.a: $(headers) $(def)/NodesHelper.cxx $(libdir)
	$(G++) -o ./$(libdir)/NodesHelper.a -c $(def)/NodesHelper.cxx

$(libdir):
	mkdir $(libdir)

$(bindir):
	mkdir $(bindir)

clean:
	\rm $(objects) && \rm $(release) temperature_*.txt

.PHONY: clean
