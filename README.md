## Compilation & Run

#### Using Makefile (recommended)

Compile only:
```
make
```

Compile and run:
```
make run
```

Compile in debug mode (with -g -G flags):
```
make debug
```

Clean compiled files:
```
make clean
```

#### Manual compilation
```
nvcc main.cu -o main && ./main
```

#### Tests
```
make test
```

```
make test ARGS="512"
```