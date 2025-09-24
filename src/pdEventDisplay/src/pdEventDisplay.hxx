#ifndef pdEventDisplay_h
#define pdEventDisplay_h

#include "pdDataClasses.hxx"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TLine.h"
#include "TMarker.h"
#include "TText.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TEllipse.h"
#include <string>
#include <vector>
#include <map>

// Structure to hold vertex information for display
struct DisplayVertex {
    double position[3];
    std::vector<AnaParticlePD*> particles;
    AnaParticlePD* parent;
    int nParticles;
    bool isTrueVertex;
    int vertexType; // 0: reconstructed, 1: true
};

class pdEventDisplay {
public:
    pdEventDisplay();
    virtual ~pdEventDisplay();

    // Initialize the event display with parameters
    bool Initialize(bool createEventDisplay = false,
                    bool saveToRootFile = false,
                    const std::string& outputDirectory = "/plots",
                    int maxEventsToDisplay = 100,
                    double eventPercentage = 100.0,
                    const std::vector<int>& requiredParticlePDGs = std::vector<int>(),
                    double vertexRadius = 30.0,
                    int minVertexDaughters = 2);

    // Main method to create event display
    void CreateEventDisplay(AnaEventB& event, int eventNumber);

    // Check if event display should be created
    bool ShouldCreateEventDisplay() const { return _CreateEventDisplay; }

    // Check if this specific event should be saved based on percentage
    bool ShouldSaveThisEvent(int eventNumber) const;

    // Check if event contains required particle types
    bool EventContainsRequiredParticles(AnaEventB& event) const;

private:
    // Event display parameters
    bool _CreateEventDisplay;
    bool _SaveToRootFile;
    std::string _OutputDirectory;
    int _MaxEventsToDisplay;
    double _EventPercentage;
    std::vector<int> _RequiredParticlePDGs;
    double _VertexRadius;
    int _MinVertexDaughters;

    // Particle type colors for display
    std::map<int, int> _particleColors;

    // Event display methods
    void InitializeParticleColors();
    void DrawEventProjections(AnaEventB& event, int eventNumber);
    void SaveCanvasToRootFile(TCanvas* canvas, int eventNumber);
    std::vector<DisplayVertex> CreateVertices(std::vector<AnaParticlePD*>& validParticles, AnaTrueParticleB* beamTrue);

    // Helper methods
    int GetParticleColor(int pdg);
    std::string GetParticleName(int pdg);
    void SetupCanvas(TCanvas* canvas);
    void DrawLegend(TLegend* legend);
    std::string ConvertProcessToString(AnaTrueParticleB::ProcessEnum process);

    // Static counter for events displayed
    static int _eventsDisplayed;
};

#endif
