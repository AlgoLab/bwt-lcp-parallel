WARNING_FLAGS := -Wall -Wno-unused-result
CFLAGS := -O2 -march=native $(WARNING_FLAGS) -DNDEBUG -D_GNU_SOURCE -std=c99 -g3
CXXFLAGS := -O2 -march=native $(WARNING_FLAGS) -DNDEBUG -D_GNU_SOURCE -std=c++11 -g3
#CXXFLAGS := $(CXXFLAGS) -UNDEBUG -O0 -g3
LDFLAGS := -lz -pthread

P := bwt_lcp decode_bwt
SUBDIRS := $(addprefix tests/, arrays supportBWT supportLists supportLCP)
TEST_INPUTS := $(addprefix tests/, paper_example.fasta bigpaper_example.fasta test.small.fasta.gz test.fasta.gz test.odd.fasta.gz test.253.fasta.gz test.290.fasta.gz)
TEST_SUMS := $(addsuffix .sums, $(TEST_INPUTS))

# Make must run serially (even if -j N > 1 is given) (GNU make extension)
.NOTPARALLEL:

.PHONY: bin
bin: $(P)

# General rule for building the executables
%: %.cpp
	$(CXX) $(CXXFLAGS) $(filter %.cpp,$^) -o $@ $(LDFLAGS)
%: %.c
	$(CC) $(CFLAGS) $(filter %.c,$^) -o $@ $(LDFLAGS)
# Additional requirements for bwt_lcp
bwt_lcp: base_types.hpp common_types.hpp encoding.hpp encoding.cpp io.hpp kseq.h

.PHONY: clean
clean:
	@echo "Cleaning..."
	rm -rf $(P) supportfile-*.tmp outfile-*.bin

.PHONY: test
test: clean bin $(TEST_INPUTS)

.PHONY: $(TEST_INPUTS)
$(TEST_INPUTS): clean bin
	@echo "Starting test on file $@ (no output)" && \
	time -v -o .btime ./bwt_lcp $@ > .log 2> .err ; \
  sed 's/^\s\+//' .btime | grep -v ': 0$$' > .time; \
	diff .log $@.log; \
	diff .err $@.err; \
	diff -y -W 120 .time $@.time; \
  rm -f .btime .time .log .err supportfile-*.tmp; \
	( test ! -s $(@).sums || md5sum -c $@.sums; ) && \
	rm -f outfile-*.bin;

.PHONY: rebuild-test
rebuild-test: clean bin $(TEST_SUMS)

.PHONY: $(TEST_SUMS)
$(TEST_SUMS): clean bin
	@echo "Rebuilding test $@ (no output)" && \
	time -v -o .time ./bwt_lcp $(basename $(@)) > $(basename $(@)).log 2> $(basename $(@)).err && \
	md5sum outfile-BWT.bin outfile-LCP.bin > $(@) || \
	true > $(@) ; \
  sed 's/^\s\+//' .time | grep -v ': 0$$' > $(basename $(@)).time ; \
  rm -f .time supportfile-*.tmp outfile-*.bin
