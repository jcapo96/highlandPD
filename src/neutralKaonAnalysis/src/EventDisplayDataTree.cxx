#include "EventDisplayDataTree.hxx"
#include "OutputManager.hxx"
#include "pdDataClasses.hxx"
#include "ToyBoxNeutralKaon.hxx"

//********************************************************************
EventDisplayDataTree::EventDisplayDataTree() : _treeInitialized(false) {
//********************************************************************
}

//********************************************************************
void EventDisplayDataTree::InitializeTree(OutputManager& output) {
//********************************************************************
    if (_treeInitialized) return; // Already initialized

    // Delegate to pdEventDisplay
    _neutralKaonEventDisplay.InitializeTree(output);

    _treeInitialized = true;
}

//********************************************************************
void EventDisplayDataTree::FillEventDisplayData(OutputManager& output, const AnaEventB& event, const ToyBoxNeutralKaon& box) {
//********************************************************************
    // Explicitly call through base class pointer to ensure virtual dispatch works
    EventDisplayBase* basePtr = &_neutralKaonEventDisplay;
    basePtr->FillTree(output, event, (void*)&box);
}
