

// Needs Sungbin's Repo
#include "SCECorr.h"

#include "Math/Vector3D.h"
using ROOT::Math::XYZVector;

SCECorr load_sce(bool is_data) {

    SCECorr* sce_corr = nullptr;
    printf("wiremod_ndhist: loading SCE TH3\n");
    sce_corr = new SCECorr(is_data);
    sce_corr->ReadHistograms();
    printf("wiremod_ndhist: loaded SCE TH3 complete\n");
    
    return sce_corr;
}



double apply_sce(SCECorr* sce_corr, XYZVector sp,

    double x,
    double y,
    double z,
    double dirx,
    double diry,
    double dirz,
    int ip

) {

    // SCE correction
    sp = sce_corr->WireToTrajectoryPosition(sp);
    double pitch_sce_uncorr = sce_corr->meas_pitch(x, y, z, dirx, diry, dirz, ip, false);
    double pitch_sce_corr = sce_corr->meas_pitch(x, y, z, dirx, diry, dirz, ip, true);
    double sce_correction = pitch_sce_uncorr / pitch_sce_corr;
    return sce_correction;

}


