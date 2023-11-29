/******************************************************************************
 *                                                                            *
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @authors: Emanuele Rambaldi <emanuelerambaldi@outlook.com>
 */

// Compile from the src folder with: g++ -g -Wall transfer_function.cpp tsv_manipulation.cpp main.cpp -o main
// For getting help on how to run the script, run from the src folder with: ./main --help
// For just reading the reference, run from the src folder with: ./main ./../data scope_data.tsv
// For both reading the reference and writing the output of the script, run from the src folder with: ./main ./../data scope_data.tsv scope_data_mod.tsv

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include "tsv_manipulation.h"
#include "transfer_function.h"

#define SAMPLING_PERIOD 0.01

int main(int argc, char *argv[]) {

  std::string extension = ".tsv";
  std::string help = "--help";
  std::string pathRd;
  std::string pathWr;
  std::vector<std::string> args;

  // Flag that is switched to true if the path to a file to write is provided
  bool wrFlag = false; 

  // Process the arguments
  if (argc == 2){
    args.push_back(argv[1]);
    if (args[0].compare(help) == 0){
      std::cout << "For just reading the reference, run from the src folder with: ./main <path-to-the-folder-of-the-tsv-file> <name-of-the-file-to-read>" << std::endl << "For both reading the reference and writing the output of the script, run from the src folder with: ./main <path-to-the-folder-of-the-tsv-file> <name-of-the-file-to-read> <name-of-the-file-to-write>" << std::endl;
      return EXIT_SUCCESS;
    }
     else {
      std::cout << "ERROR: unrecognized argument!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (argc == 3){
    args.push_back(argv[1]);
    args.push_back(argv[2]);

    std::string endingRd = args[1].substr(args[1].length()-extension.length(),extension.length());
    if (endingRd.compare(extension) == 0){
      std::cout << "The path to the .tsv file to read has been provided." << std::endl;
      // Store the path to the file that will be read
      pathRd = args[0] + "/" + args[1];

      std::cout << "Read path: " << pathRd << std::endl;
    }
    else {
      std::cout << "ERROR: the inserted name is not the name of a .tsv file!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else if (argc == 4){
    args.push_back(argv[1]);
    args.push_back(argv[2]);
    args.push_back(argv[3]);
    
    std::string endingRd = args[1].substr(args[1].length()-extension.length(),extension.length());
    std::string endingWr = args[2].substr(args[2].length()-extension.length(),extension.length());
    if ((endingRd.compare(extension) == 0) && (endingWr.compare(extension) == 0)){
      std::cout << "The paths to both the .tsv file to read and to the .tsv file to write have been provided." << std::endl;
      // Store the path to the file that will be read
      pathRd = args[0] + "/" + args[1];
      // Store the path to the file that will be written
      pathWr = args[0] + "/" + args[2];

      std::cout << "Read path: " << pathRd << std::endl;
      std::cout << "Write path: " << pathWr << std::endl;

      wrFlag = true; 
    }
    else {
      std::cout << "ERROR: the inserted name is not the name of a .tsv file!" << std::endl;
      return EXIT_FAILURE;
    }
  }
  else {
    std::cout << "ERROR: wrong number of arguments!" << std::endl;
    return EXIT_FAILURE;
  }

  // Define the vector that will contain the columns of the .tsv file
  std::vector<std::vector<double>> tsv_matrix;

  // Read the .tsv file and store the content in a matrix passed as reference
  short out = TsvManipulation::read_tsv(pathRd, tsv_matrix);
  if (out == -1) return EXIT_FAILURE;

  // Define the coefficients of both the numerator and the denominator of the transfer functions that compose the control system
  // - Discrete transfer functions
  //  - Controller
  std::vector<double> numeratorCoeffsC = {0.00045, 0.0009, 0.00045};
  std::vector<double> denominatorCoeffsC = {1.0, -1.921, 0.9231};
  //  - Reference Plant (integrator)
  double K = 1.0;
  std::vector<double> numeratorCoeffsGref = {K*SAMPLING_PERIOD/2, K*SAMPLING_PERIOD/2};
  std::vector<double> denominatorCoeffsGref = {1.0, -1.0};
  //  - Delay
  std::vector<double> numeratorCoeffsD = {0.0, 1.0};
  std::vector<double> denominatorCoeffsD = {1.0, 0.0};
  // - Continuous transfer functions
  //  - Plant
  std::vector<double> numeratorCoeffsG = {0.0, 0.0, 0.9};
  std::vector<double> denominatorCoeffsG = {0.5, 1.0, 0.0};
  //  - PI Compensator
  float Kp = 10.0;
  float Ki = 1.0;
  std::vector<double> numeratorCoeffsPI = {Kp, Ki};
  std::vector<double> denominatorCoeffsPI = {1.0, 0.0};

  // Instantiate the transfer functions that compose the control system
  // - Discrete transfer functions
  //  - Controller1
  TransferFunction dtfC1(numeratorCoeffsC, denominatorCoeffsC, SAMPLING_PERIOD, 'D');
  //  - Reference Plant (integrator)
  TransferFunction dtfGref(numeratorCoeffsGref, denominatorCoeffsGref, SAMPLING_PERIOD, 'D');
  //  - Delay1
  TransferFunction dtfD1(numeratorCoeffsD, denominatorCoeffsD, SAMPLING_PERIOD, 'D');
  //  - Controller2
  TransferFunction dtfC2(numeratorCoeffsC, denominatorCoeffsC, SAMPLING_PERIOD, 'D');
  //  - Delay2
  TransferFunction dtfD2(numeratorCoeffsD, denominatorCoeffsD, SAMPLING_PERIOD, 'D');
  // - Continuous transfer functions
  //  - Plant
  TransferFunction dtfG(numeratorCoeffsG, denominatorCoeffsG, SAMPLING_PERIOD, 'C');
  //  - PI Compensator
  TransferFunction dtfPI(numeratorCoeffsPI, denominatorCoeffsPI, SAMPLING_PERIOD, 'C');

  // Convert the continuous transfer functions to the discrete domain via the tustin bilinear transform
  // - Plant
  out = dtfG.BilinearTransform();
  if (out == -1) return EXIT_FAILURE;
  // - PI Compensator
  out = dtfPI.BilinearTransform();
  if (out == -1) return EXIT_FAILURE;
  
  // Define and Initialise inputs and outputs
  // - Reference Plant Control System
  double inputC1 = 0.0;
  double outputC1 = 0.0;
  double inputGref = 0.0;
  double outputGref = 0.0;
  std::vector<double> plant_reference; // Vector that will contain the reference plant's output signal
  double inputD1 = 0.0;
  double outputD1 = 0.0;
  // - Real Plant Control System
  double inputC2 = 0.0;
  double outputC2 = 0.0;
  double inputPI = 0.0;
  double outputPI = 0.0;
  double inputG = 0.0;
  double outputG = 0.0;
  std::vector<double> plant_output; // Vector that will contain the real plant's output signal
  double inputD2 = 0.0;
  double outputD2 = 0.0;

  // Retrieve the output of the implemented transfer function at each time step
  std::cout << "Evaluating the output of the system during " << tsv_matrix.size() << " time steps..." << std::endl;
  for (size_t i = 0; i < tsv_matrix.size(); i++) {
      // - Reference Plant Control System
      inputC1 = tsv_matrix[i][1] - outputD1;
      outputC1 = dtfC1.ComputeResponseDiscrete(inputC1);
      if (outputC1 == -1) return EXIT_FAILURE;
      inputGref = outputC1;
      outputGref = dtfGref.ComputeResponseDiscrete(inputGref);
      if (outputGref == -1) return EXIT_FAILURE;
      plant_reference.push_back(outputGref);
      inputD1 = outputGref;
      outputD1 = dtfD1.ComputeResponseDiscrete(inputD1);
      if (outputD1 == -1) return EXIT_FAILURE;
      // - Real Plant Control System
      inputC2 = tsv_matrix[i][1] - outputD2;
      outputC2 = dtfC2.ComputeResponseDiscrete(inputC2);
      if (outputC2 == -1) return EXIT_FAILURE;
      inputPI = outputD1 - outputD2;
      outputPI = dtfPI.ComputeResponseDiscrete(inputPI);
      if (outputPI == -1) return EXIT_FAILURE;
      inputG = outputC2 + outputPI;
      outputG = dtfG.ComputeResponseDiscrete(inputG);
      if (outputG == -1) return EXIT_FAILURE;
      plant_output.push_back(outputG);
      inputD2 = outputG;
      outputD2 = dtfD2.ComputeResponseDiscrete(inputD2);
      if (outputD2 == -1) return EXIT_FAILURE;
  }
  std::cout << "End of the evaluation." << std::endl;

  if(wrFlag){
    // Swap the 3-rd column of the read file with the column containing the generated plant_reference
    // Swap the 4-th column of the read file with the column containing the generated plant_output
    for (size_t i = 0; i < tsv_matrix.size(); i++) {
        tsv_matrix[i][2] = plant_reference[i];
        tsv_matrix[i][3] = plant_output[i];
    }

    // Write a .tsv file with the content of the new matrix
    out = TsvManipulation::write_tsv(tsv_matrix, pathWr);
    if (out == -1) return EXIT_FAILURE;
  }

  std::cout << "Last reference: " << plant_reference[tsv_matrix.size() - 1] << std::endl;
  std::cout << "Last output: " << plant_output[tsv_matrix.size() - 1] << std::endl;
  
  return EXIT_SUCCESS;
}
