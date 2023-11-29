Feedback Control System in C++
=================================================
# feedback-control-system_cpp

## Description
This repository contains one c++ script and two c++ libraries that reproduce the following feedback control system (instead implemented in simulink).

## Organisation
Specifically the repository is organised in 5 folders:
* `data`: folder that contains three tabular files:
	* `reference.tsv`: two-columns table containing the time steps (time: 40 seconds, sampling frequency: 100 Hz) and the corresponding input reference values
	* `reference.tsv`: four-columns table containing the time steps and the corresponding input reference values, plant reference values and plant output values (obtained from the simulink model)
	* `reference.tsv`: four-columns table containing the time steps and the corresponding input reference values, plant reference values and plant output values (obtained from the c++ script)
* `images`: folder that contains the images of this README
* `simulation`: folder that contains two MATLAB scripts and a simulink model for simulation purposes:
	* `simulink_load_tsv.m`: MATLAB script that loads the reference.tsv table into the MATLAB workspace in order to make it available for the simulink model
	* `controller.slx`: simulink model of the considered feedback control system
	* `simulink_saveplot_tsv.m`: MATLAB script that saves and plots the output data that has been saved into the MATLAB workspace by the simulink model
* src: folder that contains the main c++ script and the two c++ libraries that implement the considered control system:
	* `tsv_manipulation.cpp` and `tsv_manipulation.h`: c++ library that contains a class that in its turn is composed of two functions: one for reading tabular `.tsv` files, the other for writing tabular `.tsv` files
	* `transfer_function.cpp` and `transfer_function.h`: c++ library that contains a class that implements whatever discrete transfer function (if the object is a continuous transfer function, the class contains a method that applies the Tustin bilinear transform to turn it into the discrete domain)
	* `main.cpp`: c++ script that makes use of the two libraries to realise the considered control system
* visualisation: folder that contains a MATLAB script for visualisation purposes:
	- `plot_cpp_tsv.m`: MATLAB script that plots both the signals obtained from the reference simulink model and the ones obtained instead from the c++ scripts for comparison

## How to run
Hereafter the flow for running all the scripts of the project is presented:
* From the simulation folder run the script `simulink_load_tsv.m` on MATLAB: this will load the input reference (reference.tsv), needed for running the simulink model, into the MATLAB workspace
	<p align="center">
		<img src="" width="900" />
	</p>
* From the simulation folder run the simulink model `controller.slx`: this will produce an output (signals to be used as reference for the signals generated by the C++ script) and save the output into the MATLAB workspace
	<p align="center">
		<img src="" width="900" />
	</p>
* From the simulation folder run the script `simulink_saveplot_tsv.m`: this will both save as a `.tsv` file (`scope_data.tsv`) and plot the output of the simulink model
* Now run the c++ main script to prove that the output is the same as the one of the simulink model. To do that:
	* Compile the script from the src folder with: 
	```bash
	g++ -g -Wall transfer_function.cpp tsv_manipulation.cpp main.cpp -o main
	```
	* Run the executable that has been produced by the compilation. This could be done in 3 different ways:
		* For getting help on how to run the script, run from the src folder with: 
		```bash
		./main --help
		```
		* For just reading the reference, run from the src folder with: ./main ./../data scope_data.tsv
		```bash
		./main ./../data scope_data.tsv
		```
		* For both reading the reference and writing the output of the script, run from the src folder with: ./main ./../data scope_data.tsv scope_data_mod.tsv
		```bash
		./main ./../data scope_data.tsv scope_data_mod.tsv
		```
	* If the third modality is chosen, the output generated by the executable will be saved into the data folder as a `.tsv` file with the name `scope_data_mod.tsv`
* Finally from the visualisation folder run the script to compare the reference signals with the obtained signals: in other words, this will plot both `scope_data.tsv` and `scope_data_mod.tsv`

## Author
Emanuele Rambaldi  
E-mail: emanuelerambaldi@outlook.com 
GitHub: LaRambla20 (https://github.com/LaRambla20)
Copyright: Istituto Italiano di Tecnologia
