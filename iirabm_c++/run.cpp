#include <vector>
#include <iostream>
#include <sstream>


#include "run.h"

std::string iirabm_model_run(const std::string& params) {
    std::vector<std::string> tokens;
    std::istringstream ss(params);
    std::string token;
    while (std::getline(ss, token, ',')) {
        tokens.push_back(token);
    }
    
    float internalParam[9];
    std::string internal_param = "[";
    for (int i = 5; i < 14; ++i) {
        internalParam[i - 5] = std::stof(tokens[i]);
        if (i > 5) {
            internal_param += ", ";
        }
        internal_param += tokens[i];
    }
    internal_param += "]";

    std::cout << "Parameters: oxyHeal: " << tokens[0] << ", IS: " << tokens[1] << 
        ", NRI: " << tokens[2] << ", NIR: " << tokens[3] << ", NC: " << tokens[4] << 
        ", internalParam: " << internal_param << std::endl;

    int result = mainSimulation(std::stof(tokens[0]), std::stoi(tokens[1]),
        std::stoi(tokens[2]), std::stoi(tokens[3]), std::stoi(tokens[4]),
        internalParam);
    
    return std::to_string(result);
}