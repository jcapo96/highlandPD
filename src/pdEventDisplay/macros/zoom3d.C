// ROOT macro for 3D zooming functions
// Load this in ROOT with: .L zoom3d.C

// Function to zoom to a specific region
void Zoom3D(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax) {
    TH3F* h3d = (TH3F*)gPad->GetPrimitive("h3d");
    if (!h3d) {
        std::cout << "No 3D histogram found in current pad!" << std::endl;
        return;
    }

    // Set ranges
    h3d->GetXaxis()->SetRangeUser(xmin, xmax);
    h3d->GetYaxis()->SetRangeUser(ymin, ymax);
    h3d->GetZaxis()->SetRangeUser(zmin, zmax);

    // Update 3D view
    TView3D* view3d = (TView3D*)gPad->GetView();
    if (view3d) {
        view3d->SetRange(xmin, ymin, zmin, xmax, ymax, zmax);
    }

    // Use the canvas's built-in redraw
    gPad->Modified();
    gPad->Update();
    std::cout << "Zoomed to: X[" << xmin << "," << xmax << "] Y[" << ymin << "," << ymax << "] Z[" << zmin << "," << zmax << "]" << std::endl;
}

// Function to reset to full view
void Reset3D() {
    TH3F* h3d = (TH3F*)gPad->GetPrimitive("h3d");
    if (!h3d) {
        std::cout << "No 3D histogram found in current pad!" << std::endl;
        return;
    }

    // Reset to full detector range
    h3d->GetXaxis()->SetRangeUser(-360, 360);
    h3d->GetYaxis()->SetRangeUser(0, 700);
    h3d->GetZaxis()->SetRangeUser(0, 700);

    // Update 3D view
    TView3D* view3d = (TView3D*)gPad->GetView();
    if (view3d) {
        view3d->SetRange(-360, 0, 0, 360, 700, 700);
    }

    // Use the canvas's built-in redraw
    gPad->Modified();
    gPad->Update();
    std::cout << "Reset to full detector view" << std::endl;
}

// Function to zoom to center region
void ZoomCenter() {
    Zoom3D(-100, 100, 200, 400, 200, 400);
}

// Function to zoom to vertex region
void ZoomVertex() {
    Zoom3D(-50, 50, 300, 500, 300, 500);
}

// Function to zoom to front corner
void ZoomFront() {
    Zoom3D(-200, 200, 0, 200, 0, 200);
}

// Function to zoom to back corner
void ZoomBack() {
    Zoom3D(-200, 200, 500, 700, 500, 700);
}

// Function to show help
void ZoomHelp() {
    std::cout << "=== 3D Zoom Functions ===" << std::endl;
    std::cout << "Zoom3D(xmin, xmax, ymin, ymax, zmin, zmax) - Custom zoom" << std::endl;
    std::cout << "ZoomCenter() - Zoom to center region" << std::endl;
    std::cout << "ZoomVertex() - Zoom to vertex region" << std::endl;
    std::cout << "ZoomFront() - Zoom to front corner" << std::endl;
    std::cout << "ZoomBack() - Zoom to back corner" << std::endl;
    std::cout << "Reset3D() - Reset to full detector view" << std::endl;
    std::cout << "ZoomHelp() - Show this help" << std::endl;
    std::cout << "=========================" << std::endl;
}
