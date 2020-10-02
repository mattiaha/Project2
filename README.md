# Project2

In order to run this program, Armadillo must be installed. You must also place all the files in the same directory.
To run the program you need to run the file main.cpp, which will use every class that is included. 
When you run main.cpp the terminal will ask for the number of grid points, N. As the run-time gets significantly long when N is high, the maximum value for N allowed is N=100.
You can then choose between three functions for the rotational algorithm to solve, and you write in the number on the terminal to choose which one you want.
The program will then write results in a txt-file, both analytical and calculated. The calculated results come from our algorithm and from the armadillo library's eigenvalue solver. 
You can then use the python functions plot,plot2... etc to plot the different values produced.
plot.py lets you plot the calculated and analytical results from the txt files. You will have to write in the name of the txt file and value of N yourself.
plot2.py works the same way, but plots the relative error for the function of one electron.
plot3.py was made to show difference between rho = 5 and rho = 10 in the case of 1 electron.
plot4.py plots the ground states from the third function.
plot5.py plots the relative error between the ground states for different values of rho.

If you want to change the value of rho from 10, this can be done in the file EigenvalueSolver.cpp, in the last function QM_Matrix. 
There are two unit tests, one which makes sure the eigenvalues in the first function is within a tolerance of the analytical values.
The second test was meant to show orthogonality between eigenvectors after the algorithm was finished, but for this to work you have to change the function orthogonal in the file Test.cpp, and set the values compared to tolerance to the absolute value. 
