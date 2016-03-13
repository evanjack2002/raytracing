EXEC = raytracing
.PHONY: all
all: $(EXEC)

CC ?= gcc
CFLAGS = \
	-std=gnu99 -Wall -O0 -g
LDFLAGS = \
	-lm

ifeq ($(strip $(PROFILE)),1)
PROF_FLAGS = -pg
CFLAGS += $(PROF_FLAGS)
LDFLAGS += $(PROF_FLAGS) 
endif

ifeq ($(strip $(TEST_U)),1)
CFLAGS += -D_UNROLLING_
endif

ifeq ($(strip $(TEST_M)),1)
CFLAGS += -D_MACRO_
endif

ifeq ($(strip $(TEST_S)),1)
# -mavx -ftree-vectorize -ftree-vctorizer-verbose=1 -march=corei7-avx
CFLAGS += -D_SIMD_ -mavx
endif

OBJS := \
	objects.o \
	raytracing.o \
	main.o

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<


$(EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: use-models.h
use-models.h: models.inc Makefile
	@echo '#include "models.inc"' > use-models.h
	@egrep "^(light|sphere|rectangular) " models.inc | \
	    sed -e 's/^light /append_light/g' \
	        -e 's/light[0-9]/(\&&, \&lights);/g' \
	        -e 's/^sphere /append_sphere/g' \
	        -e 's/sphere[0-9]/(\&&, \&spheres);/g' \
	        -e 's/^rectangular /append_rectangular/g' \
	        -e 's/rectangular[0-9]/(\&&, \&rectangulars);/g' \
	        -e 's/ = {//g' >> use-models.h

gprof:
	gprof ./raytracing | ../gprof2dot/gprof2dot.py -o gr.dot
	dot -Tpng -o goutput.png gr.dot

perf:
	sudo sh -c " echo 0 > /proc/sys/kernel/kptr_restrict"
	perf record -g -- ./raytracing
	perf script | c++filt | ../gprof2dot/gprof2dot.py -f perf -o pr.dot
	dot -Tpng -o poutput.png pr.dot

rpt:
	sudo sh -c " echo 0 > /proc/sys/kernel/kptr_restrict"
	echo 1 | sudo tee /proc/sys/vm/drop_caches
#	perf stat \
#		-e cache-misses,cache-references,instructions,cycles,branches,branch-misses \
#		./$(EXEC) 1>/dev/null
	perf record \
		-e cache-misses,cache-references,instructions,cycles,branches,branch-misses \
		./$(EXEC)
	perf report

test: test.c
	$(CC) $(CFLAGS) -o $@ $<

run: $(EXEC)
	./$(EXEC)
	convert out.ppm out.png

clean:
	$(RM) $(EXEC) $(OBJS) use-models.h \
		out.ppm gmon.out *.png *.dot test
