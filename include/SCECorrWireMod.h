

// Needs Sungbin's Repo
#include "SCECorr.h"



SCECorr load_sce(bool is_data) {

    SCECorr* sce_corr = nullptr;
    if (do_sce) {
        printf("wiremod_ndhist: loading SCE TH3\n");
        sce_corr = new SCECorr(is_data);
        sce_corr->ReadHistograms();
        printf("wiremod_ndhist: loaded SCE TH3 complete\n");
    }
    return sce_corr
}



double apply_sce() {

    // SCE correction
    sp = sce_corr->WireToTrajectoryPosition(sp);
    double pitch_sce_uncorr = sce_corr->meas_pitch(x[ip][i], y[ip][i], z[ip][i], dirx[ip][i], diry[ip][i], dirz[ip][i], ip, false);
    double pitch_sce_corr = sce_corr->meas_pitch(x[ip][i], y[ip][i], z[ip][i], dirx[ip][i], diry[ip][i], dirz[ip][i], ip, true);
    dqdx_hit *= pitch_sce_uncorr / pitch_sce_corr;
    sce_correction = pitch_sce_uncorr / pitch_sce_corr;

}


