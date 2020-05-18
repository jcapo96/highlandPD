//minimal implementetation of the space charge effects in protodune
//since the most used method is voxelized_TH3, only this one is taken into account

#ifndef pdSpaceCharge_h
#define pdSpaceCharge_h

#include "TFile.h"
#include "TH3.h"

namespace pdspacecharge {
  class pdSpaceCharge {
    
  public:
    pdSpaceCharge();
    virtual ~pdSpaceCharge();
 
    bool EnableSimSpatialSCE() const;
    bool EnableSimEfieldSCE() const;
    bool EnableCalSpatialSCE() const;
    bool EnableCalEfieldSCE() const;
    bool EnableCorrSCE() const {return (EnableCalSpatialSCE()||EnableCalEfieldSCE()) ;}
    
    std::vector<double> GetPosOffsets      (std::vector<double> const& point) const;
    std::vector<double> GetEfieldOffsets   (std::vector<double> const& point) const;
    std::vector<double> GetCalPosOffsets   (std::vector<double> const& point, int const& TPCid) const;
    std::vector<double> GetCalEfieldOffsets(std::vector<double> const& point, int const& TPCid) const;
    
  private:
  protected:
         
    std::vector<double> GetOffsetsVoxel(std::vector<double> const& point, TH3F* hX, TH3F* hY, TH3F* hZ) const;
    
    std::vector<TH3F*> SCEhistograms = std::vector<TH3F*>(12); //Histograms are Dx, Dy, Dz, dEx/E0, dEy/E0, dEz/E0 (positive; repeat for negative)
    std::vector<TH3F*> CalSCEhistograms = std::vector<TH3F*>(12); 
                  
    double TransformX(double xVal) const;
    double TransformY(double yVal) const;
    double TransformZ(double zVal) const;

    bool IsInsideBoundaries(std::vector<double> const& point) const;
    bool IsTooFarFromBoundaries(std::vector<double> const& point) const;
    std::vector<double> PretendAtBoundary(std::vector<double> const& point) const;
    
    bool fEnableSimSpatialSCE;
    bool fEnableSimEfieldSCE;
    bool fEnableCalSpatialSCE;
    bool fEnableCalEfieldSCE;
    bool fEnableCorrSCE;
          
    double fEfield;
       
    std::string fInputFilename;
    std::string fCalInputFilename;

    TFile * infile;
        
    TH3F* hDx_sim_pos_orig = new TH3F();
    TH3F* hDy_sim_pos_orig = new TH3F();
    TH3F* hDz_sim_pos_orig = new TH3F();
    TH3F* hEx_sim_pos_orig = new TH3F();
    TH3F* hEy_sim_pos_orig = new TH3F();
    TH3F* hEz_sim_pos_orig = new TH3F();
    				    
    TH3F* hDx_sim_neg_orig = new TH3F();
    TH3F* hDy_sim_neg_orig = new TH3F();
    TH3F* hDz_sim_neg_orig = new TH3F();
    TH3F* hEx_sim_neg_orig = new TH3F();
    TH3F* hEy_sim_neg_orig = new TH3F();
    TH3F* hEz_sim_neg_orig = new TH3F();
    				      
    TH3F* hDx_sim_pos = new TH3F();
    TH3F* hDy_sim_pos = new TH3F();
    TH3F* hDz_sim_pos = new TH3F();
    TH3F* hEx_sim_pos = new TH3F();
    TH3F* hEy_sim_pos = new TH3F();
    TH3F* hEz_sim_pos = new TH3F();
    			       
    TH3F* hDx_sim_neg = new TH3F();
    TH3F* hDy_sim_neg = new TH3F();
    TH3F* hDz_sim_neg = new TH3F();
    TH3F* hEx_sim_neg = new TH3F();
    TH3F* hEy_sim_neg = new TH3F();
    TH3F* hEz_sim_neg = new TH3F();

    TH3F* hDx_cal_pos_orig = new TH3F();
    TH3F* hDy_cal_pos_orig = new TH3F();
    TH3F* hDz_cal_pos_orig = new TH3F();
    TH3F* hEx_cal_pos_orig = new TH3F();
    TH3F* hEy_cal_pos_orig = new TH3F();
    TH3F* hEz_cal_pos_orig = new TH3F();
    
    TH3F* hDx_cal_neg_orig = new TH3F();
    TH3F* hDy_cal_neg_orig = new TH3F();
    TH3F* hDz_cal_neg_orig = new TH3F();
    TH3F* hEx_cal_neg_orig = new TH3F();
    TH3F* hEy_cal_neg_orig = new TH3F();
    TH3F* hEz_cal_neg_orig = new TH3F();
      				     
    TH3F* hDx_cal_pos = new TH3F();
    TH3F* hDy_cal_pos = new TH3F();
    TH3F* hDz_cal_pos = new TH3F();
    TH3F* hEx_cal_pos = new TH3F();
    TH3F* hEy_cal_pos = new TH3F();
    TH3F* hEz_cal_pos = new TH3F();
    
    TH3F* hDx_cal_neg = new TH3F();
    TH3F* hDy_cal_neg = new TH3F();
    TH3F* hDz_cal_neg = new TH3F();
    TH3F* hEx_cal_neg = new TH3F();  
    TH3F* hEy_cal_neg = new TH3F();
    TH3F* hEz_cal_neg = new TH3F();
  };
} 
#endif 
