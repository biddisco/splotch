# common variables and rules for all Makefiles of RAY++

RAYROOT =       /afs/mpa/home/martin/raypp
CXX =		g++
LINKER =	g++
MAKEDEP =	g++ -E -MM -MG -I/afs/mpa/home/martin/raypp
MAKE =		gmake
CXXFLAGS =	-I/afs/mpa/home/martin/raypp -W -Wall -Wshadow -Wredundant-decls -Woverloaded-virtual -Wcast-qual -Wcast-align -Wpointer-arith -Wconversion -Wold-style-cast -fstrict-aliasing -O2 -fomit-frame-pointer -ffast-math -s  -c
LFLAGS =	
SHARELFLAGS =	-shared -o $(SHAREDLIB)
OBJECTS =	$(SOURCES:.cxx=.o)
ARCHIVE =	$(RAYROOT)/ray++.o
SHAREDLIB =	$(RAYROOT)/ray++.so
ENABLE_SHARED =	NO
ENABLE_STRIPPING = YES

.SILENT:

.PHONY: objclean depclean clean checkin compmsg ocomp compile

%.o: %.cxx
	@echo "  Compiling" $*.cxx
	$(CXX) $(CXXFLAGS) $*.cxx

.depend:
	if test -n "$(SOURCES)"; then \
		$(MAKEDEP) $(SOURCES) > .depend; \
	fi

objclean:
	rm -f $(OBJECTS)

depclean:
	rm -f .depend

clean:
	@echo Cleaning up in `pwd`
	rm -f $(OBJECTS) $(DELETABLES) .depend

compile: $(OBJECTS)

-include .depend
