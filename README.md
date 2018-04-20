# bwt_lcp

Multi-threaded construction of multistring BWT and LCP array.

## Input and Output

The input for this program should be a FASTA file on the alphabet `{A, C, G, T}`.
Some examples are provided in the directory `tests/`.  
The output of the program are the BWT and the LCP array of the collection of strings in the input file.
After the main command execution, the binary-encoded BWT is contained in the file `outfile-BWT.bin`,
while the binary LCP array is contained in the file `outfile-LCP.bin`.
Each entry of the latter file is a 8-bit unsigned integer with the only exception of `-1` encoded as `0xFF`.

## Compiling and running

After cloning the repository, issue the following command:

```
make bin
```

(Any modern C++ compiler should suffice for compiling the program.)

Now it is possible to run the program as follows:

```
./bwt_lcp <path to the FASTA file> <no of threads>
```

To output the ASCII version of the BWT use the `decode_bwt` tool as follows:

```
./decode_bwt <path to outfile-BWT.bin>
```


## Test instances

Some test instances are provided in the repository.
To run them, one can issue the command `make test`.

## Citation

Paola Bonizzoni, Gianluca Della Vedova, Yuri Pirola, Marco
Previtali and Raffaella Rizzi, _Divide and Conquer computation of the
multi-string BWT and LCP array_, to appear in CiE 2018