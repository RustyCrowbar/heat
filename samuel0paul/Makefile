headers = ./header/*.h
files = ./*cpp ./test/*.cpp ./definition/*.cxx
objects = ./lib/Nodes.a ./lib/NodesHelper.a
debug_objects = ./lib/debug_Nodes.a ./lib/debug_NodesHelper.a
Ldir = -L/usr/lib/x86_64-linux-gnu
libs = -lboost_regex
def = ./definition/
G++ = g++ -std=c++1y -Og -Wall -pthread -Wl,--no-as-needed $(Ldir)
Gdbg = g++ -std=c++1y -Og -Wall -Wshadow -ggdb -pthread -Wl,--no-as-needed $(Ldir)
libdir = lib

bindir = bin
release_target = $(bindir)/main.out
debug_target = $(bindir)/debug_main.out

all: $(release_target)

debug: $(debug_target)

gif: $(release_target)
	./$(release_target)
	./tools/make_gif.bash

$(release_target): $(files) $(objects) main.cpp $(bindir)
	$(G++) $(libs) $(objects) -o $(release_target) main.cpp

$(debug_target): $(files) $(debug_objects) main.cpp $(bindir)
	$(Gdbg) $(libs) $(objects) -o $(debug_target) main.cpp

./lib/Nodes.a: $(headers) $(def)/Nodes.cxx $(libdir)
	$(G++) -o ./$(libdir)/Nodes.a -c $(def)/Nodes.cxx

./lib/debug_Nodes.a: $(headers) $(def)/Nodes.cxx $(libdir)
	$(Gdbg) -o ./$(libdir)/Nodes.a -c $(def)/Nodes.cxx

./lib/NodesHelper.a: $(headers) $(def)/NodesHelper.cxx $(libdir)
	$(G++) -o ./$(libdir)/NodesHelper.a -c $(def)/NodesHelper.cxx

./lib/debug_NodesHelper.a: $(headers) $(def)/NodesHelper.cxx $(libdir)
	$(Gdbg) -o ./$(libdir)/NodesHelper.a -c $(def)/NodesHelper.cxx

$(libdir):
	mkdir $(libdir)

$(bindir):
	mkdir $(bindir)

clean:
	\rm $(objects) && \rm $(release) temperature_*.txt

.PHONY: clean
