pub(crate) const FM: f32 = 0.000025;
pub(crate) const PROTON_RAD: f32 = 0.8751 * FM;
pub(crate) const BOHR_RAD: f32 = 52917.721 * FM;

// In multiples of Bohr radius, the maximum distance allowed for electrons in this model
pub(crate) const MAX_RADIUS: f32 = 100.;

// Scaled radius to draw proton
pub(crate) const SCALED_PROTON_RAD: f32 = PROTON_RAD * 1000.;