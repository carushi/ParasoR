#Makefile
VAR = MED
C++11 = c++11
CC = g++
HOMEDIR:=$(shell cd $(CURDIR)/../;pwd)
CFLAGS = -O3 -std=${C++11} -Wno-unused-function -Wno-unused-variable -D$(VAR) -DGLIBCXX_DEBUG -D"HOMEDIR=\"$(HOMEDIR)/\""
LDFLAGS = -lm 
TARGET = ParasoR
OBJS = matrix.o part_func.o main.o param.o convert.o part_func_delta.o part_func_prof.o part_func_mut.o num_operator.o plot_struct.o

.PHONY: all
all: $(TARGET)
.PHONY: test
CONSTRAINT=50
CHUNK=3
test: $(TARGET)
	sh check.sh $(CONSTRAINT) $(CHUNK)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)
.cc.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<


matrix.hh: energy_const.hh
convert.hh: num_operator.hh
param.hh: pair_mat.hh convert.hh
part_func.hh: matrix.hh param.hh centroid.hh douter_concat.hh

main.o: main.cc part_func.hh part_func.o part_func_delta.o part_func_mut.o plot_struct.o part_func_prof.o
part_func.o: part_func.cc part_func.hh centroid.hh douter_concat.hh matrix.o param.o
part_func_delta.o: part_func_delta.cc part_func.hh matrix.o param.o
part_func_mut.o: part_func_mut.cc part_func.hh matrix.o param.o
part_func_prof.o: part_func_prof.cc part_func.hh matrix.o param.o
plot_struct.o: plot_struct.cc plot_struct.hh matrix.o
matrix.o: matrix.cc matrix.hh
param.o: param.cc param.hh convert.o sequence.hh
convert.o: convert.cc convert.hh num_operator.hh
num_operator.o: num_operator.cc num_operator.hh


.PHONY: clean
clean:
	rm -f *.o ParasoR *.gch *.out  *.a
.PHONY: clear
clear:
	rm -f *.o *.gch *.a
