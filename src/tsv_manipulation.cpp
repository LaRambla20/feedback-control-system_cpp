/******************************************************************************
 *                                                                            *
 * Copyright (C) 2021 Fondazione Istituto Italiano di Tecnologia (IIT)        *
 * All Rights Reserved.                                                       *
 *                                                                            *
 ******************************************************************************/

/**
 * @authors: Emanuele Rambaldi <emanuelerambaldi@outlook.com>
 */

#include "tsv_manipulation.h"


// (TsvManipulation class) Static member function that reads a .tsv file and stores its content into a matrix
int TsvManipulation::read_tsv(const std::string& path, std::vector<std::vector<double>>& matrix) {
    // Open the .tsv file
    std::ifstream inputFile(path);
    
    // Check if the file has been opened correctly
    if (!inputFile.is_open()) {
        std::cerr << "ERROR: unable to open the file." << std::endl;
        return -1;
    }

    std::string line;

    // While it is still possible to get a row from the file (getline stops when it encounters a '\n')
    while(std::getline(inputFile, line)){
        // Convert the row to a string stream so as to analyse it like an input stream
        std::istringstream ss(line);

        double value;
        std::vector<double> row;

        // Read element by element the string stream and each time store the element into the variable value (until the end of the stream is encountered)
        while (ss >> value) {
            // Add the element contained in the variable value to the considered row
            row.push_back(value);
            // Ignore the tabulation character
            ss.ignore();
        }

        // Append the considered row to the columns matrix
        matrix.push_back(row);
    }

    // Close the file
    inputFile.close();

    return 1;
}

// (TsvManipulation class) Static member function that writes a .tsv file with the content stored into a matrix
int TsvManipulation::write_tsv(const std::vector<std::vector<double>>& matrix, const std::string& path) {
    std::ostringstream tsvContent;

    for (size_t i = 0; i < matrix.size(); i++) {
        // Extract a row from the matrix
        const std::vector<double>& row = matrix[i];
        // Insert the extracted row in the stream
        for (size_t i = 0; i < row.size(); i++) {
            tsvContent << row[i];
            
            // Add a separator unless it is the last column
            if (i < row.size() - 1) {
                tsvContent << '\t';
            }
        }
        // Insert a new line after each row except for the last one
        if (i != matrix.size()-1) tsvContent << '\n'; 
    }

    // Open the file, write the stream and close the file
    std::ofstream outputFile(path);
    if (outputFile.is_open()) {
        outputFile << tsvContent.str();
        outputFile.close();
    }
    else {
        std::cerr << "ERROR: unable to open the file." << std::endl;
        return -1;
    }

    return 1;
}