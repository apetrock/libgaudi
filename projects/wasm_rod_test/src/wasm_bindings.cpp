#include <emscripten/bind.h>
#include "rod_wasm_wrapper.hpp"

using namespace emscripten;

EMSCRIPTEN_BINDINGS(rod_simulation) {
    register_vector<float>("FloatVector");
    
    class_<RodSimulationWrapper>("RodSimulation")
        .constructor<>()
        .function("step", &RodSimulationWrapper::step)
        .function("reset", &RodSimulationWrapper::reset)
        .function("getVertices", &RodSimulationWrapper::getVertices)
        .function("getNormals", &RodSimulationWrapper::getNormals)
        .function("getTangents", &RodSimulationWrapper::getTangents)
        .function("getGrowthWeights", &RodSimulationWrapper::getGrowthWeights)
        .function("getVertexCount", &RodSimulationWrapper::getVertexCount)
        .function("getTotalLength", &RodSimulationWrapper::getTotalLength)
        .function("getCurrentFrame", &RodSimulationWrapper::getCurrentFrame)
        .function("setGrowthRate", &RodSimulationWrapper::setGrowthRate)
        .function("setConstraintStrength", &RodSimulationWrapper::setConstraintStrength);
}
