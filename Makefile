all: build

build: jp.c
	clang -Ofast jp.c main.c -lgmp

bench: bench.c
	clang -Ofast bench.c -lgmp

validate: validate.c
	clang -Ofast validate.c -lgmp
