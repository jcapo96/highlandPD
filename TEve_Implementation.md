# TEve-based 3D Event Display - Implementation Complete

## Summary

The 3D event display has been **completely reimplemented** using ROOT's **TEve** (The Event Visualization Environment) framework. This provides **native mouse wheel zoom** and professional event visualization capabilities.

## ✅ What Changed

### Before (TH3F + TPolyLine3D)
- Used `TH3F` histogram as backdrop with `"ogl"` draw option
- Drew particles as `TPolyLine3D` and `TPolyMarker3D`
- **Mouse wheel zoom didn't work** on macOS
- Required custom helper macros (`zoom_3d_helper.C`)
- Manual axis range manipulation for zoom
- Limited interactivity

### After (TEve Framework)
- Uses `TEveManager` for full event visualization
- Particles drawn as `TEveLine` objects
- Vertices drawn as `TEvePointSet` objects
- **Mouse wheel zoom works natively** on all platforms
- **No helper macros needed** - everything is built-in
- Object selection, info display, and professional UI
- Better performance and rendering

## 🎯 Key Features

### Native Zoom & Rotation
- **Mouse wheel**: Zoom in/out (works immediately!)
- **Left-click + drag**: Rotate view
- **Right-click + drag**: Pan view
- **Keyboard shortcuts**: x, y, z for axis views, r to reset, etc.

### Object Selection
- **Shift+Click** on any particle to select it
- Shows particle name, type, and properties
- Highlights selected objects

### Professional UI
- Built-in camera toolbar
- Event browser sidebar
- Multiple view options (Perspective, Orthographic, etc.)
- Save screenshots directly from menu

## 📝 Modified Files

### Core Implementation
1. **`highland/src/highland2/highlandTools/src/DrawingTools.cxx`**
   - Replaced entire 3D mode implementation
   - Now uses TEve instead of TH3F/TPolyLine3D
   - Lines 3141-3264: Complete TEve implementation

2. **`highland/src/highland2/highlandTools/src/DrawingTools.cxx`** (headers)
   - Added TEve includes:
     - `#include <TEveManager.h>`
     - `#include <TEveLine.h>`
     - `#include <TEvePointSet.h>`
     - `#include <TEveViewer.h>`
     - `#include <TEveScene.h>`

### Build Configuration
3. **`highland/CMakeLists.txt`**
   - Added `Eve` component to ROOT package requirements
   - Line 42: `find_package(ROOT REQUIRED COMPONENTS Geom EG RGL Eve)`

### Documentation & Examples
4. **`DATA/test_interactive_zoom.C`**
   - Updated to reflect TEve capabilities
   - Removed obsolete zoom helper loading
   - Added instructions for built-in controls

5. **`DATA/ZOOM_INSTRUCTIONS.txt`**
   - Completely rewritten for TEve
   - Comprehensive guide to all TEve features
   - Keyboard shortcuts, toolbar usage, troubleshooting

6. **`highlandPD/TEve_Implementation.md`** (this file)
   - Implementation documentation

### Obsolete Files (kept for reference)
7. **`DATA/zoom_3d_helper.C`**
   - No longer needed with TEve
   - Kept for backward compatibility/reference

## 🚀 Usage

### Quick Start
```bash
cd /Users/jcapo/cernbox/DUNE-IFIC/Software/HIGHLAND_NEW/DATA
source ../highland/scripts/setup.sh
source ../highlandPD/scripts/setup.sh
root -l test.root
```

```cpp
.x test_interactive_zoom.C
```

### Manual Generation
```cpp
DrawingTools draw("test.root");
draw.EvtDisplay("3D", 44164030, 19, 5393);
```

### Controls
- **Zoom**: Mouse wheel up/down
- **Rotate**: Left-click + drag
- **Pan**: Right-click + drag
- **Select**: Shift+click on objects
- **Reset view**: Press 'r' key
- **Toggle sidebar**: Press 'e' key

## 🎨 Visual Elements

### Particles
- **Protons**: Blue lines
- **Pions (π+/π-)**: Red lines
- **Muons (μ+/μ-)**: Green lines
- **K0 trajectory**: Magenta dashed line

### Vertices
- **K0 Creation Vertex**: Blue star marker
- **K0 Annihilation Vertex**: Red star marker

### Labeling
- Each particle has a unique name (e.g., "Proton #5")
- Shows PDG code in element name
- K0 candidates labeled with index

## 💡 Technical Details

### TEve Initialization
```cpp
if (!gEve) {
    TEveManager::Create();
}
```

### Scene & Viewer Creation
```cpp
TEveScene* scene = gEve->SpawnNewScene(canvasName.c_str(), canvasTitle.c_str());
TEveViewer* viewer = gEve->SpawnNewViewer(canvasName.c_str());
viewer->AddScene(scene);
```

### Drawing Particles
```cpp
TEveLine* track = new TEveLine(2);
track->SetPoint(0, x1, y1, z1);
track->SetPoint(1, x2, y2, z2);
track->SetMainColor(color);
track->SetLineWidth(2);
track->SetElementName("Proton #1");
scene->AddElement(track);
```

### Camera Setup
```cpp
TGLViewer* glviewer = gEve->GetDefaultGLViewer();
glviewer->SetCurrentCamera(TGLViewer::kCameraPerspXOZ);
glviewer->CurrentCamera().RotateRad(-0.5, 0.5);
```

### Redraw
```cpp
gEve->Redraw3D(kTRUE);
```

## ✅ Benefits

1. **Native Zoom**: Mouse wheel works on all platforms (macOS, Linux, Windows)
2. **No Helper Macros**: All functionality is built-in
3. **Object Selection**: Click particles to see their properties
4. **Professional UI**: Standard HEP visualization framework
5. **Better Performance**: Optimized rendering
6. **Multiple Views**: Can have multiple TEve windows open
7. **Easy Screenshots**: Built-in save functionality
8. **Extensible**: Easy to add more visualization features

## 🔄 Backward Compatibility

The 2D mode remains unchanged and works as before. Only the 3D mode was replaced with TEve.

API remains the same:
```cpp
draw.EvtDisplay("3D", run, subrun, event);  // Now uses TEve
draw.EvtDisplay("2D", run, subrun, event);  // Still uses TCanvas
```

## 📊 Performance

TEve is generally **faster** than TH3F + TPolyLine3D for complex events because:
- Optimized OpenGL rendering
- Better object management
- Scene graph optimization
- Hardware acceleration

## 🔍 Future Enhancements

Possible additions with TEve:
- Hit-level visualization (individual detector hits)
- Shower/cluster visualization
- Detector geometry overlay
- Event animation (time evolution)
- Advanced filtering UI
- Custom color schemes
- Energy deposition visualization

## 📚 References

- [ROOT TEve Documentation](https://root.cern/doc/master/group__TEve.html)
- [TEve Tutorials](https://root.cern/doc/master/group__tutorial__eve.html)
- [Event Display Examples](https://root.cern/doc/master/alice__esd_8C.html)

## ✨ Result

**The 3D event display now has NATIVE mouse wheel zoom that works perfectly on macOS!** 🎉

No more helper macros, no more workarounds - just professional, interactive 3D visualization out of the box.

