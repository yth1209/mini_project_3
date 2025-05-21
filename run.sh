mpiCC -o project3 project3.cpp
# mpirun -np 4 ./project3 < sample/input1.txt > output/output1.txt


# gcc -o project3 project3.cpp
./project3 < sample/input1.txt > output/output1.txt
# ./project3 < sample/input2.txt > output/output2.txt
# ./project3 < sample/input3.txt > output/output3.txt
# ./project3 < sample/input4.txt > output/output4.txt
# ./project3 < sample/input5.txt > output/output5.txt


diff -bwi output/output1.txt sample/output1.txt 
# diff -bwi output/output2.txt sample/output2.txt 
# diff -bwi output/output3.txt sample/output3.txt 
# diff -bwi output/output4.txt sample/output4.txt 
# diff -bwi output/output5.txt sample/output5.txt 