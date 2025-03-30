#pragma once
#include "Sparse.h"
#include <memory>
#include <unordered_map>

class IRestriction {
public:
    virtual sparseVector operator()(int level, sparseVector& finer) = 0;
    virtual ~IRestriction() {}
};

class fullWeighting : public IRestriction {
public:
    virtual sparseVector operator()(int level, sparseVector& finer) override {
        return finer;
    }
};

class injection : public IRestriction {
public:
    virtual sparseVector operator()(int level, sparseVector& finer) override {
        return finer;
    }
};

//=============================================================================

class IInterpolation {
public:
    virtual sparseVector operator()(int level, sparseVector& coarser) = 0;
    virtual ~IInterpolation() {}
};

class linear : public IInterpolation {
public:
    virtual sparseVector operator()(int level, sparseVector& coarser) override {};
};

class quadratic : public IInterpolation {
public:
    virtual sparseVector operator()(int level, sparseVector& coarser) override {};
};

//============================================================================

class ICycle{
public:
    virtual sparseVector operator()(int level, sparseMatrix& A, sparseVector& rhs, sparseVector& initial, int nu1, int nu2,
                                    std::unique_ptr<IRestriction> restriction, std::unique_ptr<IInterpolation> interpolation) = 0;
    virtual ~ICycle() {}
};

class VCycle : public ICycle {
public:
    virtual sparseVector operator()(int level, sparseMatrix& A, sparseVector& rhs, sparseVector& initial, int nu1, int nu2,
                                    std::unique_ptr<IRestriction> restriction, std::unique_ptr<IInterpolation> interpolation) override {};
};

class FMGCycle : public ICycle {
public:
    virtual sparseVector operator()(int level, sparseMatrix& A, sparseVector& rhs, sparseVector& initial, int nu1, int nu2,
                                    std::unique_ptr<IRestriction> restriction, std::unique_ptr<IInterpolation> interpolation) override;
};

//==============================================================================

class StrategyFactory {
public:
    static std::unique_ptr<IRestriction> createRestrictionStrategy(std::string name) {
        if (name == "fullWeighting") return std::make_unique<fullWeighting>();
        else if (name == "injection") return std::make_unique<injection>();
        throw std::invalid_argument("Invalid restriction strategy name: " + name);
    }

    static std::unique_ptr<IInterpolation> createInterpolationStrategy(std::string name) {
        if (name == "linear") return std::make_unique<linear>();
        else if (name == "quadratic") return std::make_unique<quadratic>();
        throw std::invalid_argument("Invalid interpolation strategy name: " + name);
    }

    static std::unique_ptr<ICycle> createCycleStrategy(std::string name) {
        if (name == "VCycle") return std::make_unique<VCycle>();
        else if (name == "FMGCycle") return std::make_unique<FMGCycle>();
        throw std::invalid_argument("Invalid cycle strategy name: " + name);
    }
};