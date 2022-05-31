use macroquad::prelude::*;
use crate::constants::BOHR_RAD;

pub(crate) fn draw_bohr_radius(n: u32) {
    draw_sphere_wires(Vec3::ZERO, BOHR_RAD*n.pow(2) as f32, None, GRAY);
}