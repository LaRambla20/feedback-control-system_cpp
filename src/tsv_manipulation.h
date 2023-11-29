/******************************************************************************
 *                                                                            *
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @authors: Emanuele Rambaldi <emanuelerambaldi@outlook.com>
 */

#ifndef TSV_MANIPULATION_H
#define	TSV_MANIPULATION_H


#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>

// Class that contains two functions for manipulating .tsv files
class TsvManipulation{
public:
    static int read_tsv(const std::string& path, std::vector<std::vector<double>>& matrix);
    static int write_tsv(const std::vector<std::vector<double>>& matrix, const std::string& path);
};

#endif	/* TSV_MANIPULATION_H */