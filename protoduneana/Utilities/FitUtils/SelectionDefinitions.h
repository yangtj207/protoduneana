class new_interaction_topology {
 private:
  double fEndZLow, fEndZHigh;
  double fThreshold;
  bool fCexNPi0, fSignalPastFV;
 public: 
  new_interaction_topology(double endz_low, double endz_high,
                           double threshold, bool cex_nPi0, bool sig_past_fv=false)
    : fEndZLow(endz_low), fEndZHigh(endz_high), fThreshold(threshold),
      fCexNPi0(cex_nPi0), fSignalPastFV(sig_past_fv) {}

  int operator()(int pdg, double endZ,
                 std::string process, int nPi0,
                 std::vector<int> & true_daughter_pdg,
                 std::vector<double> & true_daughter_startP/*,
                 std::vector<double> & incidentEnergies*/) {

    int topology = -1;
    if (pdg == 211) {
      if (endZ < fEndZLow/*-.49375*/) {
        topology = 4;
      }
      
      else if ((endZ > fEndZHigh)/*222.10561*/ && //After FV
               ((!fSignalPastFV) || //If we don't want to consider inel. ints past APA cut
                (fSignalPastFV && process != "pi+Inelastic"))) { //If we want to consider inel past APA
        topology = 6;
      }
      else if (process == "pi+Inelastic") {
        //daughters with & without thresholds
        bool has_pion_above_threshold = false;
        for (size_t i = 0; i < true_daughter_startP.size(); ++i) {
          if (abs(true_daughter_pdg[i]) == 211 &&
              true_daughter_startP[i] > fThreshold/*.150*/) {
            has_pion_above_threshold = true;
            break;
          }
        }

        if (has_pion_above_threshold || (!fCexNPi0 && nPi0 > 1)) {
          topology = 3;
        }
        else if (nPi0 == 0) {
          topology = 1;
        }
        else if ((fCexNPi0 && nPi0 > 0) || (!fCexNPi0 && nPi0 == 1)) {
          topology = 2;
        }
        else {
          std::cout << "warning" << std::endl;
        }
      }
      else {
        topology = 7;
      }
    }
    else if (pdg == -13) {
      topology = 5;
    }
    else {
      topology = 7;
    }

    return topology;

  }
};

class inclusive_topology {
 private:
 public: 
  inclusive_topology() {}

  int operator()(int exclusive_topology) {
    if (exclusive_topology < 4) {
      return 1;
    }
    else {
      return (exclusive_topology - 2);
    }
  }
};

/*auto selection_ID = [](bool beam_is_track, bool ends_in_APA3,
                       bool no_pion_daughter,
                       bool beam_cuts, bool has_shower) {
  if (!beam_is_track) {
    return 6;
  }

  if (!beam_cuts) {
    return 5;
  }
  
  if (!ends_in_APA3) {
    return 4;
  }

  if (!no_pion_daughter) {
    return 3;
  }

  if (has_shower) {
    return 2;
  }
  else {
    return 1;
  }
};*/

class selection_ID {
 private:
  bool fDoMichel;
 public:
  selection_ID(bool do_michel) : fDoMichel(do_michel) {};

  int operator()(bool beam_is_track, bool ends_in_APA3,
                 bool no_pion_daughter,
                 bool beam_cuts, bool has_shower, bool michel_cut) {

    if (!beam_is_track) {
      return 6;
    }

    if (!beam_cuts) {
      return 5;
    }
    
    if (!ends_in_APA3) {
      return 4;
    }

    if (fDoMichel && michel_cut) {
      return 7;
    }

    if (!no_pion_daughter) {
      return 3;
    }

    if (has_shower) {
      return 2;
    }
    else {
      return 1;
    }
  }
};

class selection_ID_inclusive {
 private:
  bool fDoMichel;
 public:
  selection_ID_inclusive(bool do_michel) : fDoMichel(do_michel) {};

  int operator()(bool beam_is_track, bool ends_in_APA3,
                 bool beam_cuts, bool michel_cut) {

    if (!beam_is_track) {
      return 4;
    }

    if (!beam_cuts) {
      return 3;
    }
    
    if (!ends_in_APA3) {
      return 2;
    }

    if (fDoMichel && michel_cut) {
      return 5;
    }

    return 1;
  }
};

auto daughter_PDG_types(const std::vector<int> bt_PDGs) {
  std::vector<int> results;
  for (const int PDG : bt_PDGs) {
    if (abs(PDG) == 211) {
      results.push_back(1);
    }
    else if (abs(PDG) == 13) {
      results.push_back(2);
    }
    else if (abs(PDG) == 2212) {
      results.push_back(3);
    }
    else if (abs(PDG) == 22) {
      results.push_back(4);
    }
    else if (abs(PDG) > 2212) {
      results.push_back(5);
    }
    else if (abs(PDG) == 11) {
      results.push_back(6);
    }
    else {
      results.push_back(7);
    }
  }
  return results;
}

auto categorize_daughters = [](
    const int beam_ID,
    const std::vector<int> bt_origins, const std::vector<int> bt_IDs,
    const std::vector<int> bt_PDGs, const std::vector<int> bt_ParIDs,
    const std::vector<int> bt_ParPDGs,
    const std::vector<int> true_daughters,
    const std::vector<int> true_grand_daughters) {
  std::vector<int> results;
  for (size_t i = 0; i < bt_origins.size(); ++i) {
    if (bt_IDs[i] == beam_ID) {
      results.push_back(1);
    }
    else if (bt_origins[i] == 2) {
      results.push_back(2);
    }
    else if (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) {
      results.push_back(10); 
    }
    else if (std::find(true_daughters.begin(), true_daughters.end(), bt_IDs[i]) !=
             true_daughters.end()) {
      if (abs(bt_PDGs[i]) == 211) {
        results.push_back(3);
      }
      else if (abs(bt_PDGs[i]) == 13) {
        results.push_back(4);
      }
      else if (bt_PDGs[i] == 2212) {
        results.push_back(5);
      }
      else if (bt_PDGs[i] == 22) {
        results.push_back(6);
      }
      else if (bt_PDGs[i] > 2212) {
        results.push_back(7);
      }
      else {
        //std::cout << bt_PDGs[i] << " " << bt_ParPDGs[i] << " " << (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) << std::endl;
        results.push_back(12);
      }
    }
    else if (std::find(true_grand_daughters.begin(), 
                       true_grand_daughters.end(), bt_IDs[i]) !=
             true_grand_daughters.end()) {
      if ((bt_PDGs[i] == 22 || abs(bt_PDGs[i]) == 11) && bt_ParPDGs[i] == 111) {
        results.push_back(11);
      }
      else {
        results.push_back(8);
      }
    }
    else if (std::find(true_grand_daughters.begin(), 
                       true_grand_daughters.end(), bt_ParIDs[i]) !=
             true_grand_daughters.end()) {
        results.push_back(9);
    }
    else {
        //std::cout << bt_PDGs[i] << " " << bt_ParPDGs[i] << " " << (abs(bt_PDGs[i]) == 11 && (abs(bt_ParPDGs[i]) == 13)) << std::endl;
        results.push_back(12);
    }
  }
  return results;
};

auto leading_p_costheta = [](const double & Px,
                          const double & Py, const double & Pz,
                          const std::vector<int> & dPDGs,
                          const std::vector<double> & dPx,
                          const std::vector<double> & dPy,
                          const std::vector<double> & dPz) {
  double costheta = -999.;
  double max_p = -999.;
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  for (size_t i = 0; i < dPDGs.size(); ++i) {
    if (dPDGs[i] != 2212) continue;

    double dP = sqrt(dPx[i]*dPx[i] + dPy[i]*dPy[i] + dPz[i]*dPz[i]);
    if (dP > max_p) {
      costheta = (dPx[i]*Px + dPy[i]*Py + dPz[i]*Pz)/(dP*P);
    }
  }
  return costheta;
};

class leading_costheta {
 private:
  int fPDG;
 public:
  leading_costheta(int pdg) : fPDG(pdg) {}

  double operator()(const double & Px,
                    const double & Py, const double & Pz,
                    const std::vector<int> & dPDGs,
                    const std::vector<double> & dPx,
                    const std::vector<double> & dPy,
                    const std::vector<double> & dPz) {
    double costheta = -999.;
    double max_p = -999.;
    double P = sqrt(Px*Px + Py*Py + Pz*Pz);
    for (size_t i = 0; i < dPDGs.size(); ++i) {
      if (dPDGs[i] != fPDG) continue;

      double dP = sqrt(dPx[i]*dPx[i] + dPy[i]*dPy[i] + dPz[i]*dPz[i]);
      if (dP > max_p) {
        max_p = dP;
        costheta = (dPx[i]*Px + dPy[i]*Py + dPz[i]*Pz)/(dP*P);
      }
    }
    return costheta;
  }
};

class leading_momentum {
 private:
  int fPDG;
 public:
  leading_momentum(int pdg) : fPDG(pdg) {}

  double operator()(const std::vector<int> & dPDGs,
                    const std::vector<double> & dP) {
    double max_p = -999.;
    //double P = sqrt(Px*Px + Py*Py + Pz*Pz);
    for (size_t i = 0; i < dPDGs.size(); ++i) {
      if (dPDGs[i] != fPDG) continue;

      if (dP[i] > max_p) {
        max_p = dP[i];
      }
    }
    return max_p;
  }
};

auto backtrack_beam = [](const std::string process,
                         const bool matched,
                         const int origin, const int PDG) {
  if (process == "primary" && matched && origin == 4 && PDG == 211) {
    return 1;
  }
  else if (process == "primary" && matched && origin == 4 && PDG == -13) {
    return 2;
  }
  else if (origin == 2) {
    return 3;
  }
  else if (process == "primaryBackground") {
    return 4;
  }
  else if (process.find("Inelastic") != std::string::npos) {
    return 5;
  }
  else if (process == "Decay") {
    return 6;
  }
  else {
    return 7;
  }
};


class truncatedMean_pos {
 private:
  double fLimit;
 public:
  truncatedMean_pos(double limit)
    : fLimit(limit) {}

  std::vector<double> operator()(
      std::vector<std::vector<double>> &vecs_dEdX) {
    //size_t size = 0;
    std::vector<double> trunc_mean;
    //std::vector<double> help_vec;
    double truncate_high = 1 - fLimit;
    int i_low = 0;
    int i_high = 0;

    //sort the dEdX vecotrs in matrix
    for(auto &&vec : vecs_dEdX){
       //size = vec.size();
       std::vector<double> help_vec;


       //check dEdX vector isn't empty!
       if(vec.empty()){
          trunc_mean.push_back(-9999.);
          continue;
       }

       else{
          std::vector<double> temp_vec;
          for (double d : vec) {
            if (d < 0.) {
              continue;
            }
            temp_vec.push_back(d);
          }
          for (double d : temp_vec) {
            if (d < 0.) {
              std::cout << d << std::endl;
            }
          }
          if (temp_vec.empty()) {
            trunc_mean.push_back(-9999.);
            continue;
          }
          //Sort Vector
          sort(temp_vec.begin(), temp_vec.end());
        
          //Discard upper and lower part of signal
          //rint rounds to integer
          i_low = rint ( temp_vec.size()*fLimit);
          i_high = rint( temp_vec.size()*truncate_high);
          
          
          //if (i_high >= temp_vec.size()) std::cout << "Warning: too high" << std::endl;
          for(int i = i_low; i </*=*/ i_high; i++){
            if (temp_vec[i] < 0) {
              std::cout << "added " << temp_vec[i] << " " << i << std::endl;
            }
            help_vec.push_back(temp_vec[i]);
          };

          //Mean of help vector

          trunc_mean.push_back(accumulate(help_vec.begin(), help_vec.end(), 0.0) / help_vec.size());
          if (trunc_mean.back() < 0 && trunc_mean.back() != -9999. && trunc_mean.back() != -999.) {
            std::cout << accumulate(help_vec.begin(), help_vec.end(), 0.0) << " " << help_vec.size() << std::endl;
            std::cout << temp_vec.size() << " " << i_low << " " << i_high << std::endl;
            for (size_t i = 0; i < help_vec.size(); ++i) {
              std::cout << "\t" << help_vec[i] << std::endl;
            }
          }
       }


    }
    return trunc_mean;
  }

};


class shower_dists {
 private:
  double fTrackScoreCut;
 public:
  shower_dists(double cut)
    : fTrackScoreCut(cut) {}
  std::vector<double> operator ()(const std::vector<double> &track_score,
                                  const std::vector<double> &shower_x,
                                  const std::vector<double> &shower_y,
                                  const std::vector<double> &shower_z,
                                  double & x, double & y, double & z) {
    std::vector<double> results;
    for(size_t i = 0; i < track_score.size(); ++i){
      if ((track_score[i] < fTrackScoreCut) &&
          (track_score[i] > 0.)) {
        double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                           std::pow((shower_y[i] - y), 2) +
                           std::pow((shower_z[i] - z), 2));
        results.push_back(dist);
      }
      else {
        results.push_back(-999.);
      }
    }

    return results;
  }
};


bool shower_dist_energy_check(double shower_dist, double shower_energy) {

  double x1 = 7.5, y1 = 880.;
  double x2 = 90., y2 = 305.;

  if (shower_energy < 140. && shower_dist < 2.5) {
    return false;
  }
  else if (shower_energy < 60. && shower_dist < 5.) {
    return false;
  }
  else if (shower_energy < 20. && shower_dist < 15.) {
    return false;
  }
  else if (shower_energy > 880.) {
    return false;
  }
  else if (shower_dist > 90.) {
    return false;
  }
  //(x-x1)(y2-y1) - (y-y1)(x2-x1)
  //(x1, y1) --> (7.5, 880)
  //(x2, y2) --> (305., 90.)
  //x --> dist, y --> energy
  else if (((shower_dist - x1)*(y2 - y1) - (shower_energy - y1)*(x2 - x1)) < 0) {
    return false;
  }
  
  return true;
}

class has_shower_dist_energy {
 private:
  double fTrackScoreCut;
 public:
  has_shower_dist_energy(double cut)
    : fTrackScoreCut(cut) {}
  bool operator()(const std::vector<double> &track_score,
                  const std::vector<double> &shower_x,
                  const std::vector<double> &shower_y,
                  const std::vector<double> &shower_z,
                  const std::vector<double> &energy,
                  double & x, double & y, double & z) {
    for(size_t i = 0; i < track_score.size(); ++i){
       double dist = sqrt(std::pow((shower_x[i] - x), 2) +
                          std::pow((shower_y[i] - y), 2) +
                          std::pow((shower_z[i] - z), 2));
       
       if ((track_score[i] < fTrackScoreCut) &&
           (track_score[i] > 0.) && shower_dist_energy_check(dist, energy[i])
           /*(dist > 5. && dist < 1000.) &&
           (energy[i] > 80. && energy[i] < 1000.)*/) {
         return true;
       }
    }

    return false;
  }
};

class isBeamType {
 private:
   bool fCheckCalo;
 public:
   isBeamType(bool check_calo)
     : fCheckCalo(check_calo) {}
   bool operator()(int reco_beam_type, std::vector<double> energies) {
     if (fCheckCalo) {
       return (energies.size() > 0 && reco_beam_type == 13);
     }

     return (reco_beam_type == 13);
   }
};

class endAPA3 {
 private:
  double fEndZCut;
 public:
  endAPA3(double cut)
    : fEndZCut(cut) {}
  
  bool operator()(double reco_beam_endZ) {
    return (reco_beam_endZ < fEndZCut);
  }
};

class secondary_noPion {
 private:
  double fTrackScoreCut;
  double fChi2Cut;
  double fdEdXLow, fdEdXMed, fdEdXHigh;
 public:
  secondary_noPion(double track_score_cut, double chi2_cut,
                   double dEdX_low, double dEdX_med, double dEdX_high)
    : fTrackScoreCut(track_score_cut),
      fChi2Cut(chi2_cut),
      fdEdXLow(dEdX_low),
      fdEdXMed(dEdX_med),
      fdEdXHigh(dEdX_high) {}
  bool operator()(
      const std::vector<double> &track_score, 
      const std::vector<int> &trackID,
      const std::vector<double> &dEdX,
      const std::vector<double> &chi2,
      const std::vector<int> &ndof) {
    for( size_t i = 0; i < track_score.size(); ++i ) {
      if ((trackID[i] != -999) && (track_score[i] > fTrackScoreCut)) {
        if (dEdX[i] < fdEdXMed/*2.8*/ && dEdX[i] > fdEdXLow/*0.5*/) {
          return false;
        }
        //else if (dEdX[i] > 2.8 && dEdX[i] < 3.4) {
        else if (dEdX[i] < fdEdXHigh) {
          if (ndof[i] > 0 && chi2[i]/ndof[i] > fChi2Cut) {
            return false;
          }
        }
      }
    }

    return true;
  }
};

auto data_beam_PID = [](const std::vector<int>& pidCandidates, const int & isMC, const int & true_PDG){
  if (isMC) {
    return (true_PDG == 211 || abs(true_PDG) == 13);
  }

  auto pid_search = std::find(pidCandidates.begin(), pidCandidates.end(), 211);
  return (pid_search != pidCandidates.end());
};

class data_BI_quality {
 private:
   bool fDoTracks;
 public:
   data_BI_quality(bool do_tracks)
     : fDoTracks(do_tracks) {}
    bool operator()(int data_BI_nMomenta, int data_BI_nTracks) {

      if (fDoTracks) {
        return (data_BI_nMomenta == 1 && data_BI_nTracks == 1);
      }
      else {
        return (data_BI_nMomenta == 1);
      }
    }
};

class beam_cut_BI {
 private:
  double fXLow, fXHigh,
         fYLow, fYHigh,
         fZLow, fZHigh,
         fCosLow;
 public:
  beam_cut_BI(double x_low, double x_high,
              double y_low, double y_high,
              double z_low, double z_high,
              double cos_low)
    : fXLow(x_low), fXHigh(x_high),
      fYLow(y_low), fYHigh(y_high),
      fZLow(z_low), fZHigh(z_high),
      fCosLow(cos_low) {}
  bool operator()(double startX,
                  double startY,   double startZ,
                  double dirX,     double dirY,
                  double dirZ,     double BI_X,
                  double BI_Y,     double BI_dirX,
                  double BI_dirY,  double BI_dirZ) {

    double deltaX = startX - BI_X;
    double deltaY = startY - BI_Y;
    double cos = BI_dirX*dirX + BI_dirY*dirY +
                 BI_dirZ*dirZ;
    if( (deltaX < fXLow) || (deltaX > fXHigh) )
      return false;

    if ( (deltaY < fYLow) || (deltaY > fYHigh) )
      return false;

    if ( (startZ < fZLow) || (startZ > fZHigh) )
      return false;

    if (cos < fCosLow)
      return false;

    return true;
  }
};

class beam_cut_TPC {
 private:
  bool fDoAngle;
  //double fXCut, fYCut, fZCut, fCosLow;
  double fXYZCut, fCosCut;
  double fMeanX, fMeanY, fMeanZ;
  double fSigmaX, fSigmaY, fSigmaZ;
  double fMeanThetaX, fMeanThetaY, fMeanThetaZ;

 public:
  beam_cut_TPC(bool do_angle, double xyz_cut, double cos_cut,
               double mean_x, double mean_y, double mean_z,
               double sigma_x, double sigma_y, double sigma_z,
               double mean_theta_x, double mean_theta_y, double mean_theta_z)
   : fDoAngle(do_angle),
     fXYZCut(xyz_cut), fCosCut(cos_cut),
     fMeanX(mean_x), fMeanY(mean_y), fMeanZ(mean_z),
     fSigmaX(sigma_x), fSigmaY(sigma_y), fSigmaZ(sigma_z),
     fMeanThetaX(mean_theta_x), fMeanThetaY(mean_theta_y),
     fMeanThetaZ(mean_theta_z) {}

  bool operator()(double calo_beam_startX, double calo_beam_startY,
                  double calo_beam_startZ, double calo_beam_endX,
                  double calo_beam_endY, double calo_beam_endZ)   {

   double diffX = calo_beam_endX - calo_beam_startX;
   double diffY = calo_beam_endY - calo_beam_startY;
   double diffZ = calo_beam_endZ - calo_beam_startZ;
   double r = sqrt(diffX*diffX + diffY*diffY + diffZ*diffZ);

   double cosTrk_thetaX = diffX / r;
   double cosTrk_thetaY = diffY / r;
   double cosTrk_thetaZ = diffZ / r;

   double cosTheta = cos(fMeanThetaX*TMath::Pi()/180.)*cosTrk_thetaX +
                     cos(fMeanThetaY*TMath::Pi()/180.)*cosTrk_thetaY +
                     cos(fMeanThetaZ*TMath::Pi()/180.)*cosTrk_thetaZ;

   if (abs((calo_beam_startX - fMeanX)/fSigmaX) > fXYZCut)
      return false;

   if (abs((calo_beam_startY - fMeanY)/fSigmaY) > fXYZCut)
      return false;

   if (abs((calo_beam_startZ - fMeanZ)/fSigmaZ) > fXYZCut)
      return false;

   if (fDoAngle && (cosTheta < fCosCut))
      return false;

   return true;
   }
};

class vertex_michel_cut {
 private:
  double fCutVal;
 public:
  vertex_michel_cut(double cut_val) : fCutVal(cut_val) {}
  bool operator()(double score, int n_hits) {
    return (n_hits > 0 && (score/n_hits) > fCutVal);
  }
  
};

class fixed_interacting_energy {
  private:
   double fEnergyFix;
  public:
   fixed_interacting_energy(double fix_val = -1.) : fEnergyFix(fix_val) {}
   double operator()(const std::vector<double> & incident_energies,
                     double interacting_energy) {
     double energy = interacting_energy;
     for (size_t k = 1; k < incident_energies.size(); ++k) {
       double deltaE = ((incident_energies)[k-1] -
                        (incident_energies)[k]);
       if (deltaE > fEnergyFix && fEnergyFix > 0.) {
         energy += deltaE; 
       }
     }
     return energy;
   }
};

class beam_P_range {
  private:
    double fRangeLow, fRangeHigh;
  public:
    beam_P_range(double range_low, double range_high)
     : fRangeLow(range_low), fRangeHigh(range_high) {}
    bool operator()(double beam_inst_P) {
      return (fRangeLow < beam_inst_P*1.e3 && beam_inst_P*1.e3 < fRangeHigh);
    }
};

class beam_XY_cuts {
  private:
    double fMeanX, fMeanY, fRadius2;
  public:
    beam_XY_cuts(double x, double y, double r)
     : fMeanX(x), fMeanY(y), fRadius2(r*r) {}

    bool operator()(double beam_inst_X, double beam_inst_Y,
                    int beam_inst_nTracks) {
      if (beam_inst_nTracks != 1) return false;

      double r2 = ((beam_inst_X - fMeanX)*(beam_inst_X - fMeanX) + 
                   (beam_inst_Y - fMeanY)*(beam_inst_Y - fMeanY));
      return (r2 < fRadius2);
    }
};

class exclude_runs {
  private:
    std::vector<int> bad_runs;
  public:
    exclude_runs(std::vector<int> runs) : bad_runs(runs) {
      std::cout << "Excluded runs: " << std::endl;
      for (int & run : bad_runs) {
        std::cout << run << std::endl;
      }
    }
    bool operator()(int run) {
      return (std::find(bad_runs.begin(), bad_runs.end(), run) ==
              bad_runs.end());
    }
};

/*auto exp_coeffs(std::vector<std::vector<double>> coeffs) {
  std::vector<std::vector<double>> results;
  //w = a*exp(b*x)
  return results;
}*/

/*
class beam_inst_P_scaled {
  private:
    double fScale;
  public:
    beam_inst_P(double scale) : fScale(scale) {}
    double operator()(double beam_inst_P) {
      return fScale*beam_inst_P;
    }
};*/
