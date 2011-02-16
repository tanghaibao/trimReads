SEQAN = .

# Link against runtime library on Linux systems
OS_NAME=$(shell uname)
ifeq ($(OS_NAME),Linux)
  LDFLAGS += -lrt
endif

CPPFLAGS += -I$(SEQAN)
CPPFLAGS += -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64
CPPFLAGS += -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0
CPPFLAGS += -O3
CPPFLAGS += -W -Wall -pedantic -Wno-variadic-macros

all: $(basename $(wildcard *.cc))
test: trimReads
	$< test.fastq
clean:
	rm -f $(basename $(wildcard *.cc))
