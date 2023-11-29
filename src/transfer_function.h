/******************************************************************************
 *                                                                            *
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @authors: Emanuele Rambaldi <emanuelerambaldi@outlook.com>
 */

#ifndef TRANSFER_FUNCTION_H
#define	TRANSFER_FUNCTION_H

#include <iostream>
#include <vector>
#include <math.h>

// Class that implements a transfer function
class TransferFunction {
private:
    std::vector<double> numeratorCoeffs;
    std::vector<double> denominatorCoeffs;
    std::vector<double> inputMemory; // Vector storing the past inputs
    std::vector<double> outputMemory; // Vector storing the past outputs
    double Ts; // Sampling Period
    char dom; // Transfer Function domain

    void ShiftMemoryRegisters(double input, double output);
    
public:
    TransferFunction(const std::vector<double>& numeratorCoefficients,
                     const std::vector<double>& denominatorCoefficients,
                     double samplingPeriod,
                     char domain);
    int BilinearTransform();                 
    double ComputeResponseDiscrete(double input);
};

#endif	/* TRANSFER_FUNCTION_H */
