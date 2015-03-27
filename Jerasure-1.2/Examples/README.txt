Example Usage of piggy_back_rs_01

./piggy_back_rs_01 4 2 1 8 10
./lrc_reed_sol_01 4 2 8 2 10

./encoder <inputfile> 6 3 piggyback_rs 8 1000 1000
./decoder <filename>

#inputfile is the file to be partitioned and encoded to n files, that will be stored in ./Coding directory.
#filename input for decoder is same as input filename even if it didn't exist, this will let decoder to pick left parts of code from ./Coding directory to do reconstruction/repair/decoding
