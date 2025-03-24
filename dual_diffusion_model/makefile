#---------------------------------------------------
# get compilation flags from root-config

ROOTCFLAGS    = $(shell $(ROOTSYS)/bin/root-config --cflags)
ROOTLIBS    = $(shell $(ROOTSYS)/bin/root-config --libs) 

#---------------------------------------------------

OS            = LINUX
CXX           = g++
CXXFLAGS      = -O3 -Wall -Wno-trigraphs -fPIC
LOCALINCLUDE  = .
ROOTINCLUDE   = $(ROOTSYS)/include
INCLUDES      = -I$(LOCALINCLUDE) -I$(ROOTINCLUDE)
LD            = g++
LDFLAGS       =
SOFLAGS       = -O -shared

# installation location; please edit it to fit your environment.
INSTBASE       =$(ROOTSYS)
# FULL PATH to the MC table ; please edit it to fit your environment.
TABLE_PATH     =/home/user1/musrfit/musrfit-dynGssEALF2-v0.1.0beta/dynGssEALF2_tbl_v0.1.0.bin


# the output from the root-config script:
CXXFLAGS      += $(ROOTCFLAGS)
LDFLAGS       +=

# some definitions: headers (used to generate *Dict* stuff), sources, objects,...
OBJS =
OBJS += dynGssEALF2.o dynGssEALF2LibraryDict.o

SHLIB = libdynGssEALF2Library.so

# make the shared lib:
#
all:  $(SHLIB)

$(SHLIB): $(OBJS)
	@echo "---> Building shared library $(SHLIB) ..."
	/bin/rm -f $(SHLIB)
	$(LD) $(OBJS) $(SOFLAGS) -o $(SHLIB) $(ROOTLIBS) -lPUserFcnBase
	@echo "done"

# clean up: remove all object file (and core files)
# semicolon needed to tell make there is no source
# for this target!
#
clean:; 	@rm -f $(OBJS) *Dict* core*
	@echo "---> removing $(OBJS)"

#
$(OBJS): %.o: %.cpp
	$(CXX) $(INCLUDES) $(CXXFLAGS) -DTABLE_PATH='"$(TABLE_PATH)"' -c $<

# Generate the ROOT CINT dictionary

dynGssEALF2LibraryDict.cpp: dynGssEALF2.h dynGssEALF2LibraryLinkDef.h
	@echo "Generating dictionary $@..."
	DYLD_LIBRARY_PATH=$(ROOTSYS)/lib
	rootcint -f $@ -c -p -I$(ROOTINCLUDE) $^

install: all
	@echo "Installing shared lib: libdynGssEALF2Library.so"
ifeq ($(OS),LINUX)
	cp -pv $(SHLIB) $(INSTBASE)/lib
	cp -pv $(LOCALINCLUDE)/*.pcm $(INSTBASE)/lib
	cp -pv $(LOCALINCLUDE)/*.h $(INSTBASE)/include
endif
