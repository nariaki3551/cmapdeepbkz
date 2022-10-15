SRCS=$(wildcard *.cpp)
DEST=$(patsubst %.cpp,%,${SRCS})

all: ${DEST}

clean:
	rm ${DEST} *~

%: %.cpp *.h
	g++ -g3 -O0 -Wall -static -o $@ $< -lntl -lgmp -lm

