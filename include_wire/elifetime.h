


// Electron Lifteime Calculations




// This returns a scale factor to apply to charges

double lifetime_correction(double x, double tau) {
    static const double v_drift = 156.267;
    double out = 1.;
    if(fabs(x) > 200.) return out;

    double this_tdrift = (200. - fabs(x)) / v_drift;
    out = 1. / exp(-1. * this_tdrift / tau);

    return out;
}



