CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# the build target executable:
MODULAR_EXAMPLE = unpacker_example
TRB_EXAMPLE = trb_unpacker_example

all: $(MODULAR_EXAMPLE) $(TRB_EXAMPLE)

$(MODULAR_EXAMPLE): unpacker_example.cpp unpacker.hpp
	$(CC) $(CFLAGS) -o $@ $<

$(TRB_EXAMPLE): trb_unpacker_example.cpp trb_unpacker.hpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	$(RM) ${MODULAR_EXAMPLE} $(TRB_EXAMPLE)
