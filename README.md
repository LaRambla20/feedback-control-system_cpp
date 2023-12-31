Feedback Control System in C++
=================================================

## Description
This repository the C++ implementation of the following feedback control system (instead implemented in simulink). The system is fed with a train of two steps. This reference signal is 40 seconds long and it is sampled at 100 Hz.
<p align="center">
	<img src="https://github.com/LaRambla20/feedback-control-system_cpp/assets/91536387/305ab9d5-68f7-4030-b61c-7766c87e90bb" width="900" />
</p>

## Organisation
Specifically the repository is organised in 5 folders:
* `data`: folder that contains three tabular files:
	* `reference.tsv`: two-columns table containing the time steps and the corresponding input reference values
	* `reference.tsv`: four-columns table containing the time steps and the corresponding input reference values, plant reference values and plant output values (obtained from the simulink model)
	* `reference.tsv`: four-columns table containing the time steps and the corresponding input reference values, plant reference values and plant output values (obtained from the C++ script)
* `images`: folder that contains the images of this README
* `simulation`: folder that contains two MATLAB scripts and a simulink model for simulation purposes:
	* `simulink_load_tsv.m`: MATLAB script that loads the reference.tsv table into the MATLAB workspace in order to make it available for the simulink model
	* `controller.slx`: simulink model of the considered feedback control system
	* `simulink_saveplot_tsv.m`: MATLAB script that saves and plots the output data that has been saved into the MATLAB workspace by the simulink model
* `src`: folder that contains the main C++ script and the two C++ libraries that implement the considered control system:
	* `tsv_manipulation.cpp` and `tsv_manipulation.h`: C++ library that contains a class that in its turn is composed of two functions: one for reading tabular `.tsv` files, the other for writing tabular `.tsv` files
	* `transfer_function.cpp` and `transfer_function.h`: C++ library that contains a class that implements whatever discrete transfer function (if the object is a continuous transfer function, the class contains a method that applies the Tustin bilinear transform to convert it to the discrete domain)
	* `main.cpp`: C++ script that makes use of the two libraries to realise the considered control system
* `visualisation`: folder that contains a MATLAB script for visualisation purposes:
	- `plot_cpp_tsv.m`: MATLAB script that plots both the signals obtained from the reference simulink model and the ones obtained instead from the C++ scripts for comparison

## How to run
Hereafter the flow for running all the scripts of the project is presented:
* From the `simulation` folder run the script `simulink_load_tsv.m` on MATLAB: this will load the input reference (`reference.tsv`), needed for running the simulink model, into the MATLAB workspace
	<p align="center">
		<img src="https://github.com/LaRambla20/feedback-control-system_cpp/assets/91536387/355a63d5-2a4b-40bc-950a-227e4f26419c" width="600" />
	</p>
* From the `simulation` folder run the simulink model `controller.slx`: this will produce an output (signals to be used as reference for the signals generated by the C++ script) and save the output into the MATLAB workspace
	<p align="center">
		<img src="https://github.com/LaRambla20/feedback-control-system_cpp/assets/91536387/defdcfa5-6fa8-49ad-ab9a-f4beaa401633" width="600" />
	</p>

* From the `simulation` folder run the script `simulink_saveplot_tsv.m`: this will both save as a `.tsv` file (`scope_data.tsv`) and plot the output of the simulink model
* Now run the C++ main script to prove that the output is the same as the one of the simulink model. To do that:
	* Compile the script from the `src` folder with: 
	```bash
	g++ -g -Wall transfer_function.cpp tsv_manipulation.cpp main.cpp -o main
	```
	* Run the executable that has been produced by the compilation. This could be done in 3 different ways:
		* For getting help on how to run the script, run from the `src` folder with: 
		```bash
		./main --help
		```
		* For just reading the reference, run from the `src` folder with: ./main ./../data scope_data.tsv
		```bash
		./main ./../data scope_data.tsv
		```
		* For both reading the reference and writing the output of the script, run from the `src` folder with: ./main ./../data scope_data.tsv scope_data_mod.tsv
		```bash
		./main ./../data scope_data.tsv scope_data_mod.tsv
		```
	* If the third modality is chosen, the output generated by the executable will be saved into the `data` folder as a `.tsv` file with the name `scope_data_mod.tsv`
* Finally from the `visualisation` folder run the script to compare the signals obtained from the reference simulink model and the ones obtained instead from the C++ scripts: in other words, this will plot both `scope_data.tsv` and `scope_data_mod.tsv`
	<p align="center">
		<img src="https://github.com/LaRambla20/feedback-control-system_cpp/assets/91536387/defdcfa5-6fa8-49ad-ab9a-f4beaa401633" width="600" />
		<img src="https://github.com/LaRambla20/feedback-control-system_cpp/assets/91536387/bcd04bdd-6f0e-4b36-b151-d131b4720729" width="600" />
	</p>
 
## Author
Emanuele Rambaldi  
E-mail: emanuelerambaldi@outlook.com   
GitHub: LaRambla20 (https://github.com/LaRambla20)  
Copyright: Istituto Italiano di Tecnologia
