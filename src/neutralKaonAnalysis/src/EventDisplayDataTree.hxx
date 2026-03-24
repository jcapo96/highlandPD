#ifndef EventDisplayDataTree_h
#define EventDisplayDataTree_h

#include "neutralKaonEventDisplay.hxx"

// Forward declarations
class OutputManager;
class AnaEventB;
class ToyBoxNeutralKaon;

/// Compatibility wrapper class for event display data tree management
/// This class provides a thin wrapper around neutralKaonEventDisplay for backward compatibility
/// with existing analysis code that uses EventDisplayDataTree directly
class EventDisplayDataTree {
public:
    EventDisplayDataTree();
    virtual ~EventDisplayDataTree(){}

    /// Initialize the EventDisplayData tree (called from DefineMicroTrees)
    /// Delegates to neutralKaonEventDisplay::InitializeTree
    void InitializeTree(OutputManager& output);

    /// Fill the tree with event data
    /// Delegates to neutralKaonEventDisplay::FillTree
    void FillEventDisplayData(OutputManager& output, const AnaEventB& event, const ToyBoxNeutralKaon& box);

private:
    /// Neutral kaon analysis-specific event display implementation
    neutralKaonEventDisplay _neutralKaonEventDisplay;

    /// Flag to track if tree has been initialized
    bool _treeInitialized;
};

#endif
