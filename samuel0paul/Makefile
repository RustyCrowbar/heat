CXX_AUX = g++ -std=c++1y -Og -Wall -pthread -Wl,--no-as-needed
MKDIR=mkdir -p

objects = ./lib/Nodes.a ./lib/NodesHelper.a ./lib/config_parse.a
debug_objects = ./lib/debug_Nodes.a ./lib/debug_NodesHelper.a ./lib/config_parse.a

libs = -lboost_regex -ltbb_debug
definitions = -DTBB_USE_ASSERT -DTBB_USE_PERFORMANCE_WARNINGS

defdir = definition
libdir = lib
bindir = bin

release_target = $(bindir)/main.out
debug_target = $(bindir)/debug_main.out

all: main

main: CXX=${CXX_AUX}
main: $(bindir) $(libdir) $(release_target)

debug: CXX = ${CXX_AUX} -Wshadow -ggdb
debug: $(bindir) $(libdir) $(debug_target)

gif: main
	${RM} temperature_*.txt
	./$(release_target) -w 50 -h 50 -e 1 -p
	@./tools/make_gif.bash

$(release_target): $(objects) main.cpp
	${CXX} $(definitions) $(libs) $^ -o $@

$(debug_target): $(debug_objects) main.cpp
	${CXX} $(definitions) $(libs) $^ -o $@

$(libdir)/debug_%.a $(libdir)/%.a: $(defdir)/%.cxx
	${CXX} -c $^ -o $@

$(bindir) $(libdir):
	${MKDIR} $@

clean:
	${RM} -r $(bindir) $(libdir) temperature_*.txt

.PHONY: debug main clean
