use macroquad::math::{Vec3, vec3};
use std::f32::consts::TAU;

pub(crate) fn spherical_to_cartesian(radius: f32, polar: f32, azimuthal: f32) -> Vec3 {
    let sin_polar = polar.sin();
    let cos_polar = polar.cos();
    let sin_azimuth = azimuthal.sin();
    let cos_azimuth = azimuthal.cos();
    vec3(sin_azimuth * cos_polar, sin_polar, cos_azimuth * cos_polar) * radius
}

pub(crate) fn wrap_angle(mut angle: f32) -> f32 {
    while angle >= TAU {
        angle -= TAU;
    }
    while angle < -TAU {
        angle += TAU;
    }
    angle
}

const MAX_APPROX_EPSILON: f32 = 0.1;
// Returns max value, not max position
pub fn approx_max<F>(f: F, mut a: f32, mut b: f32, n: u32) -> f32 where F: Fn(f32) -> f32 {
    loop {
        let delta = (b - a) / n as f32;
        let mut max_pos = a;
        let mut max_val = f(a);
        for i in 1..=n {
            let pos = a + i as f32*delta;
            let val = f(pos);
            if val > max_val {
                max_pos = pos;
                max_val = val;
            }
        }
        if max_pos > a {
            a = max_pos - delta;
        }
        if max_pos < b {
            b = max_pos + delta;
        }
        if (b - a) < MAX_APPROX_EPSILON {
            return (b + a) / 2.
        }
    }
}

pub fn approx_integral<F>(f: F, a: f32, b: f32, n: u32) -> f32 where F: Fn(f32) -> f32 {
    if n % 2 == 1 {
        panic!("n must be even for Simpson's approximation");
    }
    let delta = (b - a) / n as f32;
    let mut sum = f(a) + f(b);
    for i in 0..n {
        if i % 2 == 0 {
            sum += 4. * f(a+i as f32*delta);
        } else {
            sum += 2. * f(a+i as f32*delta);
        }
    }
    delta / 3. * sum
}