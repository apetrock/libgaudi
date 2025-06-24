#pragma once

#include "gaudi/duchamp/rod_constraints_test.hpp"
#include <vector>

/**
 * Thin wrapper around the existing block_test class for WASM export
 * Does NOT modify any existing C++ code - just provides a clean interface
 */
class RodSimulationWrapper {
public:
    RodSimulationWrapper();
    ~RodSimulationWrapper() = default;

    // Simulation control
    void step(int frame);
    void reset();
    
    // Data extraction methods
    std::vector<float> getVertices() const;
    std::vector<float> getNormals() const;
    std::vector<float> getTangents() const;
    std::vector<float> getGrowthWeights() const;
    
    // Simulation state
    int getVertexCount() const;
    float getTotalLength() const;
    int getCurrentFrame() const;
    
    // Parameter control (optional - for future use)
    void setGrowthRate(float rate);
    void setConstraintStrength(float strength);

private:
    // The original simulation - kept as smart pointer, no changes needed
    gaudi::duchamp::block_test::ptr simulation_;
    int current_frame_;
    
    // Helper methods to extract data from existing simulation
    std::vector<float> extractVertexData() const;
    std::vector<float> extractNormalData() const;
};
