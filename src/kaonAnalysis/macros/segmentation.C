const int iscan = 200;
const double slength = 14;

// define the parametric line equation
void line(double t, const double *p, double &x, double &y, double &z) {
  // a parametric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1;
  x = p[0] + p[1]*t;
  y = p[2] + p[3]*t;
  z = t;
}

//function Object to be minimized
struct SumDistance2 {
  // the TGraph is a data member of the object
  TGraph2D *fGraph;
  SumDistance2(TGraph2D *g) : fGraph(g) {}
  // calculate distance line-point
  double distance2(double x,double y,double z, const double *p) {
    // distance line point is D= | (xp-x0) cross  ux |
    // where ux is direction of line and x0 is a point in the line (like t = 0)
    ROOT::Math::XYZVector xp(x,y,z);
    ROOT::Math::XYZVector x0(p[0], p[2], 0. );
    ROOT::Math::XYZVector x1(p[0] + p[1], p[2] + p[3], 1. );
    ROOT::Math::XYZVector u = (x1-x0).Unit();
    double d2 = ((xp-x0).Cross(u)).Mag2();
    return d2;
  }
  // implementation of the function to be minimized
  double operator() (const double *par) {
    assert(fGraph != 0);
    double * x = fGraph->GetX();
    double * y = fGraph->GetY();
    double * z = fGraph->GetZ();
    int npoints = fGraph->GetN();
    double sum = 0;
    for (int i  = 0; i < npoints; ++i) {
      double d = distance2(x[i],y[i],z[i],par);
      sum += d;
    }
    /*if (first) {
      std::cout << "Total Initial distance square = " << sum << std::endl;
    }
    first = false;*/
    return sum;
  }
};



void analizeTrack(int npoints, 
		  double* x, double* y, double* z){

  //draw the 3d track
  TGraph2D* tg = new TGraph2D(npoints,x,y,z);
  tg->SetMarkerStyle(20);
  tg->Draw("lp");
  //gPad->WaitPrimitive();

  //get a segment
  double vdiff[3] = {0};
  double diff     = 0;
  for(int i = 0; i < npoints; i++){
    vdiff[0] = x[i]-x[0];
    vdiff[1] = y[i]-y[0];
    vdiff[2] = z[i]-z[0];
    diff = sqrt(pow(vdiff[0],2)+pow(vdiff[1],2)+pow(vdiff[2],2));
    std::cout << diff << std::endl;
    if(diff>slength){
      std::cout << "segment" << std::endl;
      TGraph2D* tgs = new TGraph2D(i,x,y,z);
      tgs->SetMarkerStyle(20);
      tgs->SetMarkerColor(2);
      tgs->Draw("lpsame");
      //gPad->WaitPrimitive();

      // fit the graph now
      ROOT::Fit::Fitter  fitter;
      // make the functor objet
      SumDistance2 sdist(tgs);
      ROOT::Math::Functor fcn(sdist,4);
      // set the function and the initial parameter values
      double pStart[4] = {1,1,1,1};
      fitter.SetFCN(fcn,pStart);
      // set step sizes different than default ones (0.3 times parameter values)
      for (int i = 0; i < 4; ++i) fitter.Config().ParSettings(i).SetStepSize(0.01);
      bool ok = fitter.FitFCN();
      if (!ok) {
	Error("line3Dfit","Line3D Fit failed");
	return 1;
      }
      const ROOT::Fit::FitResult & result = fitter.Result();
      std::cout << "Total final distance square " << result.MinFcnValue() << std::endl;
      result.Print(std::cout);
      //gr->Draw("p0");
      // get fit parameters
      const double * parFit = result.GetParams();
      // draw the fitted line
      int n = 1000;
      double t0 = 0;
      double dt = 10;
      TPolyLine3D *l = new TPolyLine3D(n);
      for (int i = 0; i <n;++i) {
	double t = t0+ dt*i/n;
	double x,y,z;
	line(t,parFit,x,y,z);
	l->SetPoint(i,x,y,z);
      }
      l->SetLineColor(kRed);
      l->Draw();

      gPad->Update();gPad->WaitPrimitive();
      break;
    }
  }
}

void segmentation(){
  
  //the input file name
  std::string filename= "/data4/DUNE/migue/analysis/6GeV_prod4_3.root";

  //open root file
  TFile* rfile = TFile::Open(filename.c_str());
  TTree* mc    = (TTree*)rfile->Get("ana");

  //set seltrk_ndau variable
  int seltrk_ndau;
  mc->SetBranchAddress("seltrk_ndau", &seltrk_ndau);

  //two canvas
  //TCanvas* c1 = new TCanvas("c1","c1",600,600);
  //TCanvas* c2 = new TCanvas("c2","c2",1200,600);
  //c2->Divide(3,1);

  int nentries    = mc->GetEntries();
  int intentries  = 0;
  int readentries = 0;
  
  for(int i = 0; i < 1000; i += iscan){
    if(i%100==0)std::cout << i << " entries read" << std::endl;
    // is it worth looking at that bunch of entries?
    mc->Draw("seltrk_dau_truepdg","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==2 && seltrk_dau_ndau>0","",iscan,i);
    intentries = mc->GetSelectedRows();
    std::cout << intentries << " interesting entries in this bunch" << std::endl;
    if(intentries > 0){
      //if so, scan all of them
      for(int j = i; j < i+iscan; j++){
	mc->Draw("seltrk_dau_truepdg","seltrk_truepdg==211 && seltrk_dau_truepdg==321 && seltrk_dau_trueendproc==2 && seltrk_dau_ndau>0",
		 "",1,j);
	if(mc->GetSelectedRows() > 0){
	  //scan all daughters of each interesting entry
	  mc->GetEntry(j);
	  for(int k = 0; k < seltrk_ndau; k++){
	    std::stringstream ssk;
	    ssk << k;
	    mc->Draw(("seltrk_dau_hit_z["+ssk.str()+"]:seltrk_dau_hit_y["+ssk.str()+"]:seltrk_dau_hit_x["+ssk.str()+"]").c_str(),
		     ("seltrk_dau_hit_x["+ssk.str()+"]!=-999 && seltrk_truepdg==211 && seltrk_dau_truepdg["+ssk.str()+"]==321 && seltrk_dau_trueendproc["+ssk.str()+"]==2 && seltrk_dau_ndau["+ssk.str()+"]>0 && seltrk_dau_trueendmom["+ssk.str()+"]<0.4").c_str(),
		     "",1,j);
	    //is this the daughter we are looking for?
	    if(mc->GetSelectedRows() > 0){
	      //if so, analyze it
	      analizeTrack(mc->GetSelectedRows(),mc->GetV3(),mc->GetV2(),mc->GetV1());
	    } 
	  }
	  readentries++;
	  std::cout << "interesting entries read " << readentries << std::endl;
	}
	if(readentries==intentries)break;
      }
    }
    readentries = 0;
  }  
}
