SRCDIR := src
BINDIR := bin

CXXFLAGS := -std=c++20 -O3 -Wall -Wextra -march=native -I$(SRCDIR) -MMD -MP \
	$(CXXFLAGS)
LDFLAGS  := -lm

EXENAME := gagatt luminosity
EXESRCS := $(addprefix $(SRCDIR)/, $(addsuffix .cc, $(EXENAME)))
TARGETS := $(addprefix $(BINDIR)/, $(EXENAME))

LIBSRCS := $(filter-out $(EXESRCS), $(wildcard $(SRCDIR)/*.cc))
LIBOBJS := $(LIBSRCS:.cc=.o)

# GSL
CXXFLAGS += $(shell gsl-config --cflags)
LDFLAGS  += $(shell gsl-config --libs)

all: $(TARGETS)

$(BINDIR)/%: $(SRCDIR)/%.o $(LIBOBJS) | $(BINDIR)
	$(CXX) $(LDFLAGS) -o $@ $^

$(BINDIR):
	mkdir -p $@

clean:
	rm -rf $(LIBOBJS) $(BINDIR)

-include $(wildcard $(OBJDIR)/*.d)

.PHONY: all clean
