SRCDIR := src
BINDIR := bin

CXXFLAGS := -std=c++20 -O3 -Wall -Wextra -march=native -I$(SRCDIR) $(CXXFLAGS)
LDFLAGS  := -lm

EXENAME := gagatt
EXESRCS := $(addprefix $(SRCDIR)/, $(addsuffix .cc, $(EXENAME)))
TARGETS := $(addprefix $(BINDIR)/, $(EXENAME))

LIBSRCS := $(filter-out $(EXESRCS), $(wildcard $(SRCDIR)/*.cc))
LIBOBJS := $(LIBSRCS:.cc=.o)

all: $(TARGETS)

$(BINDIR)/%: $(SRCDIR)/%.o $(LIBOBJS) | $(BINDIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BINDIR):
	mkdir -p $@

clean:
	rm -rf $(LIBOBJS) $(BINDIR)

.PHONY: all clean
