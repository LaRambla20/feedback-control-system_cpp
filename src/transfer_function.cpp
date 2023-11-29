/******************************************************************************
 *                                                                            *
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @authors: Emanuele Rambaldi <emanuelerambaldi@outlook.com>
 */

#include "transfer_function.h"

// (TransferFunction class) Private method that stores the past inputs and the past outputs that are needed to evaluate the current output
void TransferFunction::ShiftMemoryRegisters(double input, double output) 
{
    inputMemory.insert(inputMemory.begin(), input); // put the current input in the first cell of the vector
    inputMemory.resize(numeratorCoeffs.size() - 1, 0.0); // resize to the size of the vector containing the numerator's coefficients -1 (thus removing any element that goes beyond that size) and fill the empty spots with 0.0

    outputMemory.insert(outputMemory.begin(), output);
    outputMemory.resize(denominatorCoeffs.size() - 1, 0.0);
}

// (TransferFunction class) Constructor
TransferFunction::TransferFunction(const std::vector<double>& numeratorCoefficients,
                                   const std::vector<double>& denominatorCoefficients,
                                   double samplingPeriod,
                                   char domain)
: numeratorCoeffs(numeratorCoefficients), denominatorCoeffs(denominatorCoefficients), Ts(samplingPeriod), dom(domain)
{
    inputMemory.resize(numeratorCoeffs.size() - 1, 0.0); // resize to the size of the vector containing the numerator's coefficients -1 and fill the empty spots with 0.0
    outputMemory.resize(denominatorCoeffs.size() - 1, 0.0);
}

// (TransferFunction class) Method that converts a continuous transfer function up to the second-order to the discrete domain via the tustin bilinear transform
int TransferFunction::BilinearTransform()
{ 
    if(dom == 'C' || dom == 'c') {
        // Define a useful constant
        float k = 2/Ts;

        // Define and Initialise the degree of both the numerator and the denominator to the number of the corresponding coefficients - 1
        short numeratorDegree = numeratorCoeffs.size() - 1;
        short denominatorDegree = denominatorCoeffs.size() - 1;
        
        // Extract the degree of both the numerator and the denominator
        size_t i = 0;
        size_t j = 0;

        for (i = 0; i < numeratorCoeffs.size(); i++){
            if(numeratorCoeffs[i] == 0) numeratorDegree--;
            else break;
        }

        for (j = 0; j < denominatorCoeffs.size(); j++){
            if(denominatorCoeffs[j] == 0) denominatorDegree--;
            else break;
        }

        // Evaluate the coefficients of the discrete transfer function based on the degrees of numerator and denominator
        if (denominatorDegree == 0){
            double c2 = denominatorCoeffs[j];

            denominatorCoeffs[j] = c2;

            if(numeratorDegree == 0){
                c2 = numeratorCoeffs[i];

                numeratorCoeffs[i] = c2;
            }
            else{
                std::cerr << "ERROR: the system is not causal: the output depends on future input values" << std::endl;
                return -1;
            }   
        }
        else if (denominatorDegree == 1){
            double c1 = denominatorCoeffs[j];
            double c2 = denominatorCoeffs[j+1];

            denominatorCoeffs[j] = c1*k + c2;
            denominatorCoeffs[j+1] = -c1*k + c2;

            if (numeratorDegree == 0 ){
                c2 = numeratorCoeffs[i];

                numeratorCoeffs[i-1] = c2;
                numeratorCoeffs[i] = c2;
            }
            else if(numeratorDegree == 1){
                c1 = numeratorCoeffs[i];
                c2 = numeratorCoeffs[i+1];

                numeratorCoeffs[i] = c1*k + c2;
                numeratorCoeffs[i+1] = -c1*k + c2;
            }
            else{
                std::cerr << "ERROR: the system is not causal: the output depends on future input values" << std::endl;
                return -1;
            }
        }
        else if (denominatorDegree == 2){
            double c0 = denominatorCoeffs[j];
            double c1 = denominatorCoeffs[j+1];
            double c2 = denominatorCoeffs[j+2];

            denominatorCoeffs[j] = c0*powf(k,2) + c1*k + c2;
            denominatorCoeffs[j+1] = -2*c0*powf(k,2) + 2*c2;
            denominatorCoeffs[j+2] = c0*powf(k,2) - c1*k + c2;

            if (numeratorDegree == 0 ){
                c2 = numeratorCoeffs[i];

                numeratorCoeffs[i-2] = c2;
                numeratorCoeffs[i-1] = 2*c2;
                numeratorCoeffs[i] = c2;
            }
            else if(numeratorDegree == 1){
                c1 = numeratorCoeffs[i];
                c2 = numeratorCoeffs[i+1];

                numeratorCoeffs[i-1] = c1*k + c2;
                numeratorCoeffs[i] = 2*c2;
                numeratorCoeffs[i+1] = -c1*k + c2;
            }
            else if(numeratorDegree == 2){
                c0 = numeratorCoeffs[i];
                c1 = numeratorCoeffs[i+1];
                c2 = numeratorCoeffs[i+2];

                numeratorCoeffs[i] = c0*powf(k,2) + c1*k + c2;
                numeratorCoeffs[i+1] = -2*c0*powf(k,2) + 2*c2;
                numeratorCoeffs[i+2] = c0*powf(k,2) - c1*k + c2;
            }
            else{
                std::cerr << "ERROR: the system is not causal: the output depends on future input values" << std::endl;
                return -1;
            }
        }
        else{
            std::cerr << "ERROR: the function converts transfer functions that are up to the second order!" << std::endl;
            return -1;
        }
        
        // Keep track of the fact that the function has been converted before returning
        dom = 'D';
        return 0;
    }
    else {
        std::cerr << "ERROR: the input transfer function is not continuous!" << std::endl;
        return -1;
    }
}

// (TransferFunction class) Method that evaluates the output of the discrete transfer function at the current time step by implementing the corresponding difference equation
double TransferFunction::ComputeResponseDiscrete(double input)
{   
    if(dom == 'D' || dom == 'd') {
        double output = input * numeratorCoeffs[0];

        for (size_t i = 1; i < numeratorCoeffs.size(); i++) {
            output += numeratorCoeffs[i] * inputMemory[i-1];
        }

        for (size_t i = 1; i < denominatorCoeffs.size(); i++) {
            output -= denominatorCoeffs[i] * outputMemory[i-1];
        }

        output /= denominatorCoeffs[0];

        // Update the memory registers that contain the past inputs and the past outputs
        ShiftMemoryRegisters(input, output);

        return output;
    }
    else {
        std::cerr << "ERROR: the input transfer function is not discrete!" << std::endl;
        return -1.0;
    }
}