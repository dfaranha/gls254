CC=clang
ifeq "$(CC)" "gcc"
    FLAGS = -march=armv8.4-a+crypto -O3
else ifeq "$(CC)" "clang"
    FLAGS = -O3
endif

runsandbox:
	$(CC) -o run_sandbox $(FLAGS) sandboxmain.c common/*.c common/linear_pass.S
	./run_sandbox

runtests:
	$(CC) -o run_tests $(FLAGS) testmain.c common/*.c tests/*.c benchmark/benchmark_tool.c common/linear_pass.S
	./run_tests

runbench:
	$(CC) -o run_bench $(FLAGS) benchmain.c common/*.c benchmark/*.c common/linear_pass.S
	./run_bench