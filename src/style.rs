use macroquad::prelude::*;
use macroquad::ui::{Skin, root_ui};
use crate::quantum::Wavefunction;
use num_complex::Complex32;
use crate::util::spherical_to_cartesian;
use std::fmt::{Display, Formatter, Debug};
use std::f32::consts::{PI, TAU};

pub(crate) fn regular_text_format() -> egui::TextFormat {
    egui::TextFormat {
        font_id: egui::FontId {
            size: 25.,
            ..Default::default()
        },
        ..Default::default()
    }
}
pub(crate) fn subscript_text_format() -> egui::TextFormat {
    egui::TextFormat {
        font_id: egui::FontId {
            size: 15.,
            ..Default::default()
        },
        valign: egui::Align::BOTTOM,
        ..Default::default()
    }
}
pub(crate) fn superscript_text_format() -> egui::TextFormat {
    egui::TextFormat {
        font_id: egui::FontId {
            size: 15.,
            ..Default::default()
        },
        valign: egui::Align::TOP,
        ..Default::default()
    }
}

pub(crate) fn default_skin() -> Skin {
    let label_style = root_ui()
        .style_builder()
        .text_color(WHITE)
        .font_size(30)
        .build();

    let button_style = root_ui()
        .style_builder()
        .text_color(BLACK)
        .font_size(20)
        .margin(RectOffset::new(10., 10., 10., 10.))
        .color(LIGHTGRAY)
        .color_hovered(LIGHTGRAY)
        .color_clicked(GRAY)
        .build();

    Skin {
        label_style,
        button_style,
        ..root_ui().default_skin()
    }
}

const DIST_COLOR_SCALE: [Color; 3] = [
    ORANGE, MAGENTA, BLUE
];
const PROB_COLOR_SCALE: [Color; 5] = [
    BLACK, DARKPURPLE, RED, ORANGE, YELLOW
];

#[derive(Eq, PartialEq)]
pub(crate) enum ColorScheme {
    Charge,
    Phase,
    Distance,
    Probability
}
impl ColorScheme {
    pub fn get_proton_color(&self) -> Color {
        match self {
            ColorScheme::Charge => RED,
            ColorScheme::Phase => GREEN,
            ColorScheme::Distance => LIGHTGRAY,
            ColorScheme::Probability => LIGHTGRAY,
            _ => LIGHTGRAY
        }
    }

    pub fn get_electron_color(&self, wf: &Wavefunction, max_prob: f32, r: f32, max_rad: f32, polar: f32, azimuth: f32) -> Color {
        match self {
            ColorScheme::Charge => BLUE,
            ColorScheme::Phase => if wf.total_wave(r, polar, azimuth).arg().is_sign_positive() {RED} else {BLUE},
            ColorScheme::Distance => Self::color_scale(&DIST_COLOR_SCALE, r / max_rad),
            ColorScheme::Probability => Self::color_scale(&PROB_COLOR_SCALE, wf.total_pdf(r, polar, azimuth) / max_prob)
        }
    }

    fn color_scale(scale: &[Color], intensity: f32) -> Color {
        let gradient = intensity.clamp(0., 1.) * (scale.len() - 1) as f32;
        let l = gradient % 1.;
        if l == 0. {
            scale[gradient as usize]
        } else {
            let c1 = scale[gradient.floor() as usize];
            let c2 = scale[gradient.ceil() as usize];
            Color {
                r: c1.r + l * (c2.r - c1.r),
                g: c1.g + l * (c2.g - c1.g),
                b: c1.b + l * (c2.b - c1.b),
                a: c1.a + l * (c2.a - c1.a),
            }
        }
    }
}
impl Display for ColorScheme {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(match self {
            ColorScheme::Charge => "Charge Color",
            ColorScheme::Phase => "Phase Color",
            ColorScheme::Distance => "Distance Gradient",
            ColorScheme::Probability => "Probability Heatmap",
        }, f)
    }
}