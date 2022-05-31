use crate::sampling::{rejection_sample};
use std::f32::consts::{FRAC_PI_2, PI, TAU};
use macroquad::prelude::get_internal_gl;
use scilib::math::polynomial::{Laguerre, Legendre};
use scilib::quantum::spherical_harmonics;
use scilib::math::basic::factorial;
use num_complex::Complex32;
use std::cell::RefCell;
use crate::util::approx_max;
use crate::constants::MAX_RADIUS;
use crate::style::{regular_text_format, subscript_text_format};

const ORBITAL_LETTER: [char; 7] = ['s','p','d','f','g','h','i'];

pub struct Electron {
    pub r: f32,
    pub polar: f32,
    pub azimuth: f32,
    pub total_prob: f32
}

pub struct Wavefunction {
    n: u32,
    l: u32,
    m: i32,
    abs_m: u32,
    r_norm: f32,
    laguerre: Laguerre,
    //sph_norm: f32,
    polar_norm: f32,
    azimuth_norm: f32,
    legendre: Legendre,
    r_max: f32,
    polar_max: f32,
    azimuth_max: f32
}
impl Wavefunction {
    pub fn new(n: u32, l: u32, m: i32) -> Self {
        if l >= n {
            panic!("Azimuthal quantum number (l) must be less than principal quantum number (n)")
        }
        let abs_m = m.abs() as u32;
        if abs_m > l {
            panic!("Magnetic quantum number must be between -l and l (azimuthal quantum number)")
        }
        let nf = n as f32;
        let mut wf = Self {
            n,
            l,
            m,
            abs_m,
            r_norm: ((2. / nf).powi(3) * factorial((n - l - 1) as usize) as f32 / (2. * nf * factorial((n + l) as usize) as f32)).sqrt(),
            laguerre: Laguerre::new((n-l-1) as usize, (2*l+1) as i32),
            //r_sampler_state: RefCell::new(MetropolisSamplerState::new(1.)),
            //sph_norm: (-1i32).pow(abs_m) as f32 * (((2*l+1) as f32 * factorial((l-abs_m) as usize) as f32) / (4. * PI * factorial((l+abs_m) as usize) as f32)).sqrt(),
            polar_norm: (1. / PI).sqrt(),
            azimuth_norm: (((2*l+1) as f32 * factorial((l-abs_m) as usize) as f32) / (4. * factorial((l+abs_m) as usize) as f32)).sqrt(),
            legendre: Legendre::new(l as usize, abs_m as i32),
            r_max: 0.,
            polar_max: 0.,
            azimuth_max: 0.
        };
        wf.r_max = approx_max(|r| wf.r_pdf(r), 0., MAX_RADIUS, 100);
        wf.polar_max = approx_max(|p| wf.polar_pdf(p), 0., TAU, 100);
        wf.azimuth_max = approx_max(|a| wf.azimuth_pdf(a), 0., TAU, 100);
        wf
    }

    pub fn layout_orbital_name(&self, layout: &mut egui::text::LayoutJob) {
        layout.append(&format!("{}{}", self.n, ORBITAL_LETTER[self.l as usize]), 0., regular_text_format());
        layout.append(&self.m.to_string(), 0., subscript_text_format());
    }

    pub fn get_n(&self) -> u32 {
        self.n
    }
    pub fn get_l(&self) -> u32 {
        self.l
    }
    pub fn get_m(&self) -> i32 {
        self.m
    }

    // r in multiples of Bohr radius
    pub fn r_wave(&self, r: f32) -> f32 {
        let rho = 2. * r / self.n as f32;
        self.r_norm * (-0.5 * rho).exp() * rho.powi(self.l as i32) * self.laguerre.compute(rho as f64) as f32
    }

    pub fn r_pdf(&self, r: f32) -> f32 {
        let psi = self.r_wave(r);
        (r*psi).powi(2)
    }

    pub fn polar_wave(&self, polar: f32) -> Complex32 {
        self.polar_norm * Complex32::new(0., self.m as f32 * polar).exp()
    }

    pub fn polar_pdf(&self, polar: f32) -> f32 {
        let cpx = self.polar_wave(polar);
        (cpx.re.powi(2) + cpx.im.powi(2)) as f32
    }

    /*pub fn polar_wave_real(&self, polar: f32) -> f32 {
        self.polar_norm * (self.m as f32 * polar).cos()
    }

    pub fn polar_pdf_real(&self, polar: f32) -> f32 {
        self.polar_wave_real(polar).powi(2)
    }*/

    pub fn azimuth_wave(&self, azimuth: f32) -> f32 {
        self.azimuth_norm * self.legendre.compute(azimuth.cos() as f64) as f32
    }

    pub fn azimuth_pdf(&self, azimuth: f32) -> f32 {
        self.azimuth_wave(azimuth).powi(2)
    }

    /*pub fn sph_wave(&self, polar: f32, azimuthal: f32) -> Complex32 {
        self.sph_norm * self.legendre.compute(azimuthal.cos() as f64) as f32 * Complex32::new(0., self.m as f32 * polar).exp()
    }

    pub fn sph_wave_real(&self, polar: f32, azimuthal: f32) -> f32 {
        self.sph_norm * self.legendre.compute(azimuthal.cos() as f64) as f32 * (self.m as f32 * polar).cos()
    }

    pub fn sph_pdf(&self, polar: f32, azimuthal: f32) -> f32 {
        let cpx = self.sph_wave(polar, azimuthal);
        (cpx.re.powi(2) + cpx.im.powi(2)) as f32
        //self.sph_wave_real(polar, azimuthal).powi(2)
    }*/

    pub fn total_wave(&self, r: f32, polar: f32, azimuthal: f32) -> Complex32 {
        self.r_wave(r) * self.polar_wave(polar) * self.azimuth_wave(azimuthal)
    }

    /*pub fn total_wave_real(&self, r: f32, polar: f32, azimuthal: f32) -> f32 {
        self.r_wave(r) * self.polar_wave_real(polar) * self.azimuth_wave(azimuthal)
    }*/

    pub fn total_pdf(&self, r: f32, polar: f32, azimuthal: f32) -> f32 {
        self.r_pdf(r) * self.polar_pdf(polar) * self.azimuth_pdf(azimuthal)
    }

    pub fn measure(&self) -> Electron {
        let r = rejection_sample(|r| self.r_pdf(r), self.r_max);//metropolis_sample(&self.r_sampler_state, |r| self.r_pdf(r));

        let polar = rejection_sample(|p| self.azimuth_pdf(p), self.polar_max);
        let azimuth = rejection_sample(|a| self.azimuth_pdf(a), self.azimuth_max);

        Electron {
            r,
            polar,
            azimuth,
            total_prob: self.total_pdf(r, polar, azimuth)
        }
    }
}

