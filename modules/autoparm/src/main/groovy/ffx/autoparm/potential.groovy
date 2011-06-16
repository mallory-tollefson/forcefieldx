//POTENTIAL
//Supply a choice (1 - 4)
//1. Get QM Potential from a Gaussian CUBE File
//2. Calculate the Model Potential for a System
//3. Compare a Model Potential to a Target Grid
//4. Fit Electrostatic Parameters to Target Grid
//if choice == 1, then supply the cube filename
//if choice == 2 || choice == 3, then supply the xyz filename
//if choice == 4, then supply the xyz filename and the eps


//Get XYZ File
Double eps = null;
Integer choice = Integer.parseInt(args[0]);
String filename = args[1];
if(choice == 4){
	eps = Double.parseDouble(args[2]);
}
potential(choice, filename, eps);
