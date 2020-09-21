#include <string>

extern "C" int mainSimulation(float oxyHeal, int infectSpread, int numRecurInj, 
    int numInfectRepeat, int numCytokines, float* internalParameterization);


std::string iirabm_model_run(const std::string& params);