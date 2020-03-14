ifndef TOFPET3DDIR
$(error TOFPET3DDIR is not set)
endif

OBJDIR = $(TOFPET3DDIR)/obj
SRCDIR = $(TOFPET3DDIR)/cc
LIBDIR = $(TOFPET3DDIR)/lib

SOURCES = $(wildcard $(SRCDIR)/*.cc)
OBJECTS = $(subst $(SRCDIR),$(OBJDIR),$(SOURCES:.cc=.o))
TARGETS = $(OBJECTS)
TARGETS += dirs

CXX = g++

# Compilation flags
CFLAGS = -fpic -c

# Optimization flags
CFLAGS += -O2

# Debug flags
CFLAGS += -pg

# Linking flags
#LDFLAGS = -L/Users/jrenner/IFIC/software/root/lib -lGeom -lm -lGui -L/opt/local/lib/gcc48 -lCore \

all:	$(TARGETS)
	@echo Creating library libmlem...
	@$(CXX) -shared -o $(LIBDIR)/libmlem.so $(OBJECTS)
	@echo Finished.

dirs : 
	@if [ ! -d $(LIBDIR) ] ; then mkdir -p $(LIBDIR); \
	    echo "   >>>> Created directory $(LIBDIR)"; \
	else echo " $(LIBDIR) already exists"; fi
	@if [ ! -d $(OBJDIR) ]; then mkdir -p $(OBJDIR); \
	    echo "   >>>> Create directory $(OBJDIR)"; \
	else echo " $(OBJDIR) already exists"; fi
        
clean:
	@echo Removing object files...
	@$(RM) $(OBJDIR)/*.o
	@echo Removing libraries...
	@$(RM) $(LIBDIR)/*.so

$(OBJDIR)/mlem.o: \
	$(SRCDIR)/mlem.cc
	@echo $@
	@$(CXX) $(CFLAGS) $< -o $@
