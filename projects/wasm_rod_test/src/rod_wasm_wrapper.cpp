#include "rod_wasm_wrapper.hpp"

RodSimulationWrapper::RodSimulationWrapper() 
    : current_frame_(0) {
    // Create the simulation exactly as the original code does
    simulation_ = gaudi::duchamp::block_test::create();
}

void RodSimulationWrapper::step(int frame) {
    if (simulation_) {
        simulation_->step(frame);
        current_frame_ = frame;
    }
}

void RodSimulationWrapper::reset() {
    // Recreate the simulation to reset state
    current_frame_ = 0;
    simulation_ = gaudi::duchamp::block_test::create();
}

std::vector<float> RodSimulationWrapper::getVertices() const {
    return extractVertexData();
}

std::vector<float> RodSimulationWrapper::getNormals() const {
    return extractNormalData();
}

std::vector<float> RodSimulationWrapper::getTangents() const {
    if (!simulation_ || !simulation_->__R) {
        return std::vector<float>();
    }
    
    // Extract tangent vectors from the rod
    std::vector<gaudi::vec3> tangents = simulation_->__R->N2c();
    std::vector<float> result;
    result.reserve(tangents.size() * 3);
    
    for (const auto& t : tangents) {
        result.push_back(static_cast<float>(t.x()));
        result.push_back(static_cast<float>(t.y()));
        result.push_back(static_cast<float>(t.z()));
    }
    
    return result;
}

std::vector<float> RodSimulationWrapper::getGrowthWeights() const {
    if (!simulation_ || !simulation_->__R) {
        return std::vector<float>();
    }
    
    // Get growth weights - this calls the existing method without modification
    std::vector<gaudi::real> weights = simulation_->compute_growth_weights(current_frame_);
    std::vector<float> result;
    result.reserve(weights.size());
    
    for (gaudi::real w : weights) {
        result.push_back(static_cast<float>(w));
    }
    
    return result;
}

int RodSimulationWrapper::getVertexCount() const {
    if (!simulation_ || !simulation_->__R) {
        return 0;
    }
    return static_cast<int>(simulation_->__R->__x.size());
}

float RodSimulationWrapper::getTotalLength() const {
    if (!simulation_ || !simulation_->__R) {
        return 0.0f;
    }
    return static_cast<float>(simulation_->__R->lavg() * simulation_->__R->__x.size());
}

int RodSimulationWrapper::getCurrentFrame() const {
    return current_frame_;
}

void RodSimulationWrapper::setGrowthRate(float rate) {
    // Future implementation - could modify simulation parameters
    // For now, this is a placeholder
}

void RodSimulationWrapper::setConstraintStrength(float strength) {
    // Future implementation - could modify simulation parameters
    // For now, this is a placeholder
}

std::vector<float> RodSimulationWrapper::extractVertexData() const {
    if (!simulation_ || !simulation_->__R) {
        return std::vector<float>();
    }
    
    // Extract vertex positions from the existing rod data structure
    const std::vector<gaudi::vec3>& vertices = simulation_->__R->__x;
    std::vector<float> result;
    result.reserve(vertices.size() * 3);
    
    for (const auto& v : vertices) {
        result.push_back(static_cast<float>(v.x()));
        result.push_back(static_cast<float>(v.y()));
        result.push_back(static_cast<float>(v.z()));
    }
    
    return result;
}

std::vector<float> RodSimulationWrapper::extractNormalData() const {
    if (!simulation_ || !simulation_->__R) {
        return std::vector<float>();
    }
    
    // Extract normal vectors from the rod
    std::vector<gaudi::vec3> normals = simulation_->__R->N0c();
    std::vector<float> result;
    result.reserve(normals.size() * 3);
    
    for (const auto& n : normals) {
        result.push_back(static_cast<float>(n.x()));
        result.push_back(static_cast<float>(n.y()));
        result.push_back(static_cast<float>(n.z()));
    }
    
    return result;
}
