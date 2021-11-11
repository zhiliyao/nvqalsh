SRCS=random.cc pri_queue.cc util.cc nvm-common.cc mempool.cc \
	tree.cc lbtree.cc qalsh.cc qalsh_tree.cc ann.cc \
	main.cc 
	
OBJS=${SRCS:.cc=.o}
LIB=-lpmem

CXX=g++ -std=c++11
CPPFLAGS=-w -O3
# CPPFLAGS=-w -O0 -g

.PHONY: clean

all: ${OBJS}
	${CXX} ${CPPFLAGS} -o qalsh ${OBJS} ${LIB}

random.o: random.h

pri_queue.o: pri_queue.h

util.o: util.h

nvm-common.o: nvm-common.h

mempool.o: mempool.h

tree.o: tree.h

lbtree.o: lbtree.h

qalsh.o: qalsh.h

qalsh_tree.o: qalsh_tree.h

ann.o: ann.h

main.o:

clean:
	-rm ${OBJS}
