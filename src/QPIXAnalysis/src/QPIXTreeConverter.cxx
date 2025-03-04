#include "Parameters.hxx"
#include "QPIXTreeConverter.hxx"

//********************************************************************
QPIXTreeConverter::QPIXTreeConverter(const std::string& name):InputConverter(name){
//********************************************************************

  //constructor
  _spill = NULL;
  _filename = "";
}

//********************************************************************
bool QPIXTreeConverter::Initialize(){
//********************************************************************

  AddChain(_treeName);
  fChain = GetChain(_treeName);
  
  // Set object pointer
  InitializeVariables();

  // Set branch addresses and branch pointers
  if (!fChain) return false;
  fCurrent = -1;
  SetBranchAddresses();
  
  return true;
}

//********************************************************************
QPIXTreeConverter::~QPIXTreeConverter(){
//********************************************************************
  
  if(fChain)delete fChain->GetCurrentFile();
}

//****************************************************************************
bool QPIXTreeConverter::AddFileToTChain(const std::string& inputString){
//****************************************************************************

  std::cout << "QPIXTreeConverter::AddFileToTChain(). Adding file: " << inputString << std::endl;

  _filename = inputString;

  // Chain only the directories we are interested in
  if(fChain) fChain->AddFile(inputString.c_str());
  
  // Sets the software version for this file
  return header().SetSoftwareVersion("");
}


//*****************************************************************************
Int_t QPIXTreeConverter::ReadEntry(Long64_t& entry) {
//*****************************************************************************
  
  Int_t entry_temp = fChain->GetEntry(entry);

  return entry_temp;
}

//*****************************************************************************
Int_t QPIXTreeConverter::GetSpill(Long64_t& entry, AnaSpillC*& spill){
//*****************************************************************************

  Int_t entry_temp = ReadEntry(entry);

  if(entry_temp > 0){
    
    // Create an instance of the Spill
    spill = MakeSpill();
    
    // Cast it to AnaSpill
    _spill = static_cast<AnaSpill*>(spill);
    
    // Fill the EventModel
    FillInfo(_spill);
  }
  else{
    std::cout << "Failed in reading entry " << entry << std::endl;
  }

  entry++;

  
  return entry_temp;
}

//*****************************************************************************
void QPIXTreeConverter::FillInfo(AnaSpill* spill){
//*****************************************************************************

  spill->EventInfo   = MakeEventInfo();
  spill->DataQuality = MakeDataQuality();
  spill->Beam        = MakeBeam();
  spill->Bunches.push_back(MakeBunch());

  // general event info
  FillEventInfo(static_cast<AnaEventInfoQPIX*>(spill->EventInfo));
  
  // True information
  FillTrueInfo(spill);
}

//*****************************************************************************
void QPIXTreeConverter::FillEventInfo(AnaEventInfoQPIX* info){
//*****************************************************************************

  info->Run          = _qpix_run;
  info->Event        = _qpix_event;
  info->IsBackground = _is_background;
  info->HasPotasium  = _has_potasium;
  info->HasGamma     = _has_gamma;
  for(int iplane = 0; iplane < 7; iplane++){
    info->DetectedPhotons[iplane] = _detected_photons->at(iplane);
    for(int ihit = 0; ihit < (int)_wf[iplane]->size(); ihit++)
      info->PDSwf[iplane].push_back(_wf[iplane]->at(ihit));
  }

  GetSolarChain(info);
  GetNuReaction(info);
}

//*****************************************************************************
void QPIXTreeConverter::GetSolarChain(AnaEventInfoQPIX* info){
//*****************************************************************************
  
  for(int i = AnaEventInfoQPIX::SolarChainEnum::unknown_chain; 
      i != AnaEventInfoQPIX::SolarChainEnum::dummy_chain; i++){
    if(_filename.find(info->ChainToString(i)) != std::string::npos){
      info->SolarChain = AnaEventInfoQPIX::SolarChainEnum(i);
      break;
    }
  }
}

//*****************************************************************************
void QPIXTreeConverter::GetNuReaction(AnaEventInfoQPIX* info){
//*****************************************************************************
  
  for(int i = AnaEventInfoQPIX::SolarNuReactionEnum::unknown_reaction; 
      i != AnaEventInfoQPIX::SolarNuReactionEnum::dummy_reaction; i++){
    if(_filename.find(info->NuReactionToString(i)) != std::string::npos){
      info->NuReaction = AnaEventInfoQPIX::SolarNuReactionEnum(i);
      break;
    }
  }
}

//*****************************************************************************
void QPIXTreeConverter::FillTrueInfo(AnaSpill* spill){
//*****************************************************************************

  // Clear the true particles vector
  spill->TrueParticles.clear();
  
  // loop over particles and add them to the vector of true particles
  for(int ipart = 0; ipart < _number_particles; ipart++){
    AnaTrueParticlePD* truePart = MakeTrueParticle();
    truePart->ID           = _particle_track_id->at(ipart);
    truePart->ParentID     = _particle_parent_track_id->at(ipart);
    truePart->PDG          = _particle_pdg_code->at(ipart);
    truePart->Momentum     = _particle_initial_energy->at(ipart);
    //truePart->ProcessStart = truePart->ConvertProcess(_particle_process_key->at(ipart));
    truePart->Position[0]  = _particle_initial_x->at(ipart);
    truePart->Position[1]  = _particle_initial_y->at(ipart);
    truePart->Position[2]  = _particle_initial_z->at(ipart);
    truePart->Position[3]  = _particle_initial_t->at(ipart);
    double norm = sqrt(pow(_particle_initial_px->at(ipart),2)
		       +pow(_particle_initial_py->at(ipart),2)
		       +pow(_particle_initial_pz->at(ipart),2));
    truePart->Direction[0] = _particle_initial_px->at(ipart)/sqrt(norm);
    truePart->Direction[1] = _particle_initial_py->at(ipart)/sqrt(norm);
    truePart->Direction[2] = _particle_initial_pz->at(ipart)/sqrt(norm);  
    spill->TrueParticles.push_back(truePart);    
  }

  //fill the true vertex depending on signal/background
  if(_is_background){
    if(_particle_initial_energy->empty())return; //safety
    spill->TrueVertices.clear();
    AnaTrueVertex* trueVertex = MakeTrueVertex();
    trueVertex->NuEnergy    = _particle_initial_energy->at(0);
    trueVertex->Position[0] = _particle_initial_x->at(0);
    trueVertex->Position[1] = _particle_initial_y->at(0);
    trueVertex->Position[2] = _particle_initial_z->at(0);
    trueVertex->Position[3] = _particle_initial_t->at(0);
    spill->TrueVertices.push_back(trueVertex);
    return;
  }
  else{
    if(_generator_initial_particle_energy->empty())return; //safety
    spill->TrueVertices.clear();
    AnaTrueVertex* trueVertex = MakeTrueVertex();
    //loop over initial generator particles and get neutrino energy
    for(int ipart = 0; ipart < (int)_generator_initial_particle_energy->size(); ipart++){
      if(_generator_initial_particle_pdg_code->at(ipart) == 12){
	trueVertex->NuEnergy    = _generator_initial_particle_energy->at(ipart);
	trueVertex->Position[0] = _generator_initial_particle_x->at(ipart);
	trueVertex->Position[1] = _generator_initial_particle_y->at(ipart);
	trueVertex->Position[2] = _generator_initial_particle_z->at(ipart);
	trueVertex->Position[3] = _generator_initial_particle_t->at(ipart);
      }
      if(_generator_final_particle_pdg_code->at(ipart) == 11){
	trueVertex->LeptonMom = _generator_final_particle_energy->at(ipart);
      }
    }
    spill->TrueVertices.push_back(trueVertex);  
  }
}

//*****************************************************************************
void QPIXTreeConverter::InitializeVariables(){
//*****************************************************************************

  _qpix_run = 0;
  _qpix_event = 0;
  _number_particles = 0;
  _is_background = 0;
  _has_potasium = 0;
  _has_gamma = 0;
  _particle_track_id = 0;
  _particle_parent_track_id = 0;
  _particle_pdg_code = 0;
  _particle_initial_energy = 0;
  _particle_mass = 0;
  _particle_charge = 0;
  _particle_process_key = 0;
  _particle_initial_x = 0;
  _particle_initial_y = 0;
  _particle_initial_z = 0;
  _particle_initial_t = 0;
  _particle_initial_px = 0;
  _particle_initial_py = 0;
  _particle_initial_pz = 0;
  _particle_number_daughters = 0;
  _hit_track_id = 0;
  _hit_start_x = 0;
  _hit_start_y = 0;
  _hit_start_z = 0;
  _hit_start_t = 0;
  _hit_end_x = 0;
  _hit_end_y = 0;
  _hit_end_z = 0;
  _hit_end_t = 0;
  _hit_length = 0;
  _hit_energy_deposit = 0;
  _hit_process_key = 0;
  _generator_initial_particle_energy = 0;
  _generator_initial_particle_pdg_code = 0;
  _generator_initial_particle_mass = 0;
  _generator_initial_particle_charge = 0;
  _generator_initial_particle_x = 0;
  _generator_initial_particle_y = 0;
  _generator_initial_particle_z = 0;
  _generator_initial_particle_t = 0;
  _generator_initial_particle_px = 0;
  _generator_initial_particle_py = 0;
  _generator_initial_particle_pz = 0;
  _generator_final_particle_energy = 0;
  _generator_final_particle_pdg_code = 0;
  _generator_final_particle_mass = 0;
  _generator_final_particle_charge = 0;
  _generator_final_particle_x = 0;
  _generator_final_particle_y = 0;
  _generator_final_particle_z = 0;
  _generator_final_particle_t = 0;
  _generator_final_particle_px = 0;
  _generator_final_particle_py = 0;
  _generator_final_particle_pz = 0;  
  for(int i = 0; i < 7; i++)_wf[i] = 0;
  _detected_photons = 0;
}

//*****************************************************************************
void QPIXTreeConverter::SetBranchAddresses(){
//***************************************************************************** 

  fChain->SetBranchAddress("run"                                , &_qpix_run                           , &b_qpix_run);
  fChain->SetBranchAddress("event"                              , &_qpix_event			       , &b_qpix_event);
  fChain->SetBranchAddress("number_particles"                   , &_number_particles		       , &b_number_particles);
  fChain->SetBranchAddress("is_background"                      , &_is_background		       , &b_is_background);
  fChain->SetBranchAddress("has_potasium"                       , &_has_potasium		       , &b_has_potasium);
  fChain->SetBranchAddress("has_gamma"                          , &_has_gamma	         	       , &b_has_gamma);
  fChain->SetBranchAddress("particle_track_id"                  , &_particle_track_id		       , &b_particle_track_id);
  fChain->SetBranchAddress("particle_parent_track_id"           , &_particle_parent_track_id	       , &b_particle_parent_track_id);
  fChain->SetBranchAddress("particle_pdg_code"                  , &_particle_pdg_code		       , &b_particle_pdg_code);
  fChain->SetBranchAddress("particle_initial_energy"            , &_particle_initial_energy	       , &b_particle_initial_energy);
  fChain->SetBranchAddress("particle_mass"                      , &_particle_mass		       , &b_particle_mass);
  fChain->SetBranchAddress("particle_charge"                    , &_particle_charge		       , &b_particle_charge);
  fChain->SetBranchAddress("particle_process_key"               , &_particle_process_key	       , &b_particle_process_key);
  fChain->SetBranchAddress("particle_initial_x"                 , &_particle_initial_x		       , &b_particle_initial_x);
  fChain->SetBranchAddress("particle_initial_y"                 , &_particle_initial_y		       , &b_particle_initial_y);
  fChain->SetBranchAddress("particle_initial_z"                 , &_particle_initial_z		       , &b_particle_initial_z);
  fChain->SetBranchAddress("particle_initial_t"                 , &_particle_initial_t		       , &b_particle_initial_t);
  fChain->SetBranchAddress("particle_initial_px"                , &_particle_initial_px		       , &b_particle_initial_px);
  fChain->SetBranchAddress("particle_initial_py"                , &_particle_initial_py		       , &b_particle_initial_py);
  fChain->SetBranchAddress("particle_initial_pz"                , &_particle_initial_pz		       , &b_particle_initial_pz);
  fChain->SetBranchAddress("particle_number_daughters"          , &_particle_number_daughters	       , &b_particle_number_daughters);
  fChain->SetBranchAddress("hit_track_id"                       , &_hit_track_id		       , &b_hit_track_id);
  fChain->SetBranchAddress("hit_start_x"                        , &_hit_start_x			       , &b_hit_start_x);
  fChain->SetBranchAddress("hit_start_y"                        , &_hit_start_y			       , &b_hit_start_y);
  fChain->SetBranchAddress("hit_start_z"                        , &_hit_start_z			       , &b_hit_start_z);
  fChain->SetBranchAddress("hit_start_t"                        , &_hit_start_t			       , &b_hit_start_t);
  fChain->SetBranchAddress("hit_end_x"                          , &_hit_end_x			       , &b_hit_end_x);
  fChain->SetBranchAddress("hit_end_y"                          , &_hit_end_y			       , &b_hit_end_y);
  fChain->SetBranchAddress("hit_end_z"                          , &_hit_end_z			       , &b_hit_end_z);
  fChain->SetBranchAddress("hit_end_t"                          , &_hit_end_t			       , &b_hit_end_t);
  fChain->SetBranchAddress("hit_length"                         , &_hit_length			       , &b_hit_length);
  fChain->SetBranchAddress("hit_energy_deposit"                 , &_hit_energy_deposit		       , &b_hit_energy_deposit);
  fChain->SetBranchAddress("hit_process_key"                    , &_hit_process_key		       , &b_hit_process_key);
  fChain->SetBranchAddress("generator_initial_particle_energy"  , &_generator_initial_particle_energy  , &b_generator_initial_particle_energy);
  fChain->SetBranchAddress("generator_initial_particle_pdg_code", &_generator_initial_particle_pdg_code, &b_generator_initial_particle_pdg_code);
  fChain->SetBranchAddress("generator_initial_particle_mass"    , &_generator_initial_particle_mass    , &b_generator_initial_particle_mass);
  fChain->SetBranchAddress("generator_initial_particle_charge"  , &_generator_initial_particle_charge  , &b_generator_initial_particle_charge);
  fChain->SetBranchAddress("generator_initial_particle_x"       , &_generator_initial_particle_x       , &b_generator_initial_particle_x);
  fChain->SetBranchAddress("generator_initial_particle_y"       , &_generator_initial_particle_y       , &b_generator_initial_particle_y);
  fChain->SetBranchAddress("generator_initial_particle_z"       , &_generator_initial_particle_z       , &b_generator_initial_particle_z);
  fChain->SetBranchAddress("generator_initial_particle_t"       , &_generator_initial_particle_t       , &b_generator_initial_particle_t);
  fChain->SetBranchAddress("generator_initial_particle_px"      , &_generator_initial_particle_px      , &b_generator_initial_particle_px);
  fChain->SetBranchAddress("generator_initial_particle_py"      , &_generator_initial_particle_py      , &b_generator_initial_particle_py);
  fChain->SetBranchAddress("generator_initial_particle_pz"      , &_generator_initial_particle_pz      , &b_generator_initial_particle_pz);
  fChain->SetBranchAddress("generator_final_particle_energy"    , &_generator_final_particle_energy    , &b_generator_final_particle_energy);
  fChain->SetBranchAddress("generator_final_particle_pdg_code"  , &_generator_final_particle_pdg_code  , &b_generator_final_particle_pdg_code);
  fChain->SetBranchAddress("generator_final_particle_mass"      , &_generator_final_particle_mass      , &b_generator_final_particle_mass);
  fChain->SetBranchAddress("generator_final_particle_charge"    , &_generator_final_particle_charge    , &b_generator_final_particle_charge);
  fChain->SetBranchAddress("generator_final_particle_x"         , &_generator_final_particle_x	       , &b_generator_final_particle_x);
  fChain->SetBranchAddress("generator_final_particle_y"         , &_generator_final_particle_y	       , &b_generator_final_particle_y);
  fChain->SetBranchAddress("generator_final_particle_z"         , &_generator_final_particle_z	       , &b_generator_final_particle_z);
  fChain->SetBranchAddress("generator_final_particle_t"         , &_generator_final_particle_t	       , &b_generator_final_particle_t);
  fChain->SetBranchAddress("generator_final_particle_px"        , &_generator_final_particle_px	       , &b_generator_final_particle_px);
  fChain->SetBranchAddress("generator_final_particle_py"        , &_generator_final_particle_py	       , &b_generator_final_particle_py);
  fChain->SetBranchAddress("generator_final_particle_pz"        , &_generator_final_particle_pz        , &b_generator_final_particle_pz);
  fChain->SetBranchAddress("wf_total"                           , &_wf[0]                              , &b_wf[0]);
  fChain->SetBranchAddress("wf_plane_0"                         , &_wf[1]                              , &b_wf[1]);
  fChain->SetBranchAddress("wf_plane_1"                         , &_wf[2]                              , &b_wf[2]);
  fChain->SetBranchAddress("wf_plane_2"                         , &_wf[3]                              , &b_wf[3]);
  fChain->SetBranchAddress("wf_plane_3"                         , &_wf[4]                              , &b_wf[4]);
  fChain->SetBranchAddress("wf_plane_4"                         , &_wf[5]                              , &b_wf[5]);
  fChain->SetBranchAddress("wf_plane_5"                         , &_wf[6]                              , &b_wf[6]);
  fChain->SetBranchAddress("detected_photons"                   , &_detected_photons                   , &b_detected_photons);
}
