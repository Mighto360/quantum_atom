mod constants;
mod util;
mod sampling;
mod style;
mod bohr;
mod quantum;

use crate::util::{spherical_to_cartesian, wrap_angle, approx_integral};
use crate::constants::*;
use crate::style::{default_skin, ColorScheme, regular_text_format, superscript_text_format};
use crate::bohr::draw_bohr_radius;
use crate::quantum::{Electron, Wavefunction};
use std::f32::consts::{FRAC_PI_2, PI};
use macroquad::prelude::*;
use std::ops::{Add, Mul};
use scilib::quantum::radial_wavefunction;
use scilib::constant::A_0;
use std::borrow::BorrowMut;
use egui::Pos2;
use egui::panel::TopBottomSide;
use macroquad::rand::srand;

const ANGULAR_SPEED: f32 = 1.;
const MIN_ZOOM: f32 = 3.;
const MAX_ZOOM: f32 = 175.;
const MAX_POLAR: f32 = FRAC_PI_2-0.0001;

#[macroquad::main("QuantumAtom")]
async fn main() {
    // TODO: Random seed
    srand(6032171360124809224u64);
    //let default_skin = default_skin();
    //root_ui().push_skin(&default_skin);

    let mut num_particles = 1;
    let draw_pdfs = false;
    let mut show_bohr_radius = false;
    let mut show_crosshair = false;
    let mut show_debug = false;
    let mut color_scheme = ColorScheme::Phase;

    let global_up = vec3(0., 1., 0.);
    let target = vec3(0., 0., 0.);

    let mut wf = Wavefunction::new(1, 0, 0);
    let mut electrons = Vec::<Electron>::new();
    let mut max_prob = 0.;
    let mut max_rad = 0.;

    let mut sample = "Hello";

    let mut radius = MIN_ZOOM;
    let mut polar_angle = 0.;
    let mut azimuthal_angle = 0.;
    let mut last_mouse_pos = mouse_position();
    let mut egui_using_pointer = false;
    loop {
        // Input and movement
        let delta = get_frame_time();
        let mouse_pos = mouse_position();
        if !egui_using_pointer {
            match mouse_wheel() {
                (_x, y) if y != 0.0 => {
                    // Normalize mouse wheel values is browser (chromium: 53, firefox: 3)
                    #[cfg(target_arch = "wasm32")]
                        let y = if y < 0.0 {
                        -1.0
                    } else if y > 0.0 {
                        1.0
                    } else {
                        0.0
                    };
                    radius = (radius * 1.1f32.powf(-y)).clamp(MIN_ZOOM, MAX_ZOOM);
                }
                _ => (),
            }

            if is_mouse_button_down(MouseButton::Left) {
                azimuthal_angle = wrap_angle(azimuthal_angle - (mouse_pos.0 - last_mouse_pos.0) * delta * ANGULAR_SPEED);
                polar_angle = wrap_angle(polar_angle + (mouse_pos.1 - last_mouse_pos.1) * delta * ANGULAR_SPEED).clamp(-MAX_POLAR, MAX_POLAR);
            }
        }
        last_mouse_pos = mouse_pos;

        // Scene
        clear_background(BLACK);

        let camera_pos = spherical_to_cartesian(radius, polar_angle, azimuthal_angle);
        let forward = (target - camera_pos).normalize();
        let right = forward.cross(global_up).normalize();
        let up = right.cross(forward).normalize();
        let camera = Camera3D {
            position: camera_pos,
            up: global_up,
            target,
            ..Default::default()
        };
        let camera_mat = camera.matrix();
        set_camera(&camera);

        //draw_grid(20, 1., BLACK, GRAY);

        // Nucleus
        /*draw_line_3d(vec3(-0.25, 0., 0.), vec3(0.25, 0., 0.), WHITE);
        draw_line_3d(vec3(0., 0., -0.25), vec3(0., 0., 0.25), WHITE);*/
        let proton_color = color_scheme.get_proton_color();
        /*if !show_crosshair {
            let crosshair_rad = 0.25;
            draw_line_3d(-crosshair_rad*right, crosshair_rad*right, proton_color);
            draw_line_3d(-crosshair_rad*up, crosshair_rad*up, proton_color);
        }*/
        draw_sphere(Vec3::ZERO, SCALED_PROTON_RAD, None, proton_color);
        if show_bohr_radius {
            draw_bohr_radius(wf.get_n());
        }

        //let mut i = 0;
        let quad_gl = unsafe {
            get_internal_gl().quad_gl
        };
        for e in &electrons {
            let world_pos = spherical_to_cartesian(e.r, e.polar, e.azimuth) * BOHR_RAD;
            let sides = particle_lod(camera_pos, world_pos, forward, 90.);
            let color = color_scheme.get_electron_color(&wf, max_prob, e.r, max_rad, e.polar, e.azimuth);
            draw_poly_billboard(quad_gl, up, right, world_pos, sides, 0.003*radius, 0., color);
        }

        if draw_pdfs {
            for a in 0..20 {
                let azimuth = a as f32 * PI / 10.;
                for i in 0..100 {
                    let r = i as f32 * 0.25;
                    let world_pos = spherical_to_cartesian(r, 0., azimuth) * BOHR_RAD
                        + vec3(0., 10.*wf.r_pdf(r), 0.);
                    draw_poly_billboard(quad_gl, up, right, world_pos, 10, 0.003*radius, 0., GREEN);
                }
            }
        }

        // GUI
        set_default_camera();
        let center_x = screen_width() / 2.;
        let center_y = screen_height() / 2.;
        if show_crosshair {
            let crosshair_rad = 15.;
            draw_line(center_x-crosshair_rad, center_y, center_x+crosshair_rad, center_y, 2., proton_color);
            draw_line(center_x, center_y-crosshair_rad, center_x, center_y+crosshair_rad, 2., proton_color);
        }
        let proton_scale = 1. - (radius - MIN_ZOOM) / (12. - MIN_ZOOM);
        if proton_scale > 0. {
            draw_text_ex("P+", center_x-10.*proton_scale, center_y-15.*proton_scale, TextParams {
                font_size: 15,
                font_scale: proton_scale,
                color: proton_color,
                ..Default::default()
            });
        }
        if show_debug {
            let debug_color = LIGHTGRAY;
            draw_text(&format!("Camera ({}, {:.2}°, {:.2}°)", radius, polar_angle.to_degrees(), azimuthal_angle.to_degrees()), 10., 60., 11., debug_color);
            draw_text(&format!("{} FPS ({:.0} ms)", get_fps(), delta*1000.), 10., 80., 11., debug_color);
        }

        egui_macroquad::ui(|ctx| {
            egui::TopBottomPanel::top("Top Panel").show(ctx, |ui| {
                ui.horizontal(|ui| {
                    // Left
                    /*let mut element_popup_layout = egui::text::LayoutJob::default();
                    element_popup_layout.append("H", 0., regular_text_format());
                    element_popup_layout.append("1", 0., superscript_text_format());
                    let element_popup_response = ui.add_sized(egui::Vec2::new(45., 30.), egui::Button::new(element_popup_layout));
                    let element_popup_id = ui.make_persistent_id("Element Popup Button");
                    if element_popup_response.clicked() {
                        ui.memory().toggle_popup(element_popup_id);
                    }
                    egui::popup_below_widget(ui, element_popup_id, &element_popup_response, |ui| {
                        ui.set_min_width(100.);
                    });*/
                    ui.vertical(|ui| {
                        ui.label("View");
                        ui.horizontal(|ui| {
                            ui.menu_button("Go to Location", |ui| {
                                if ui.button("Nucleus").clicked() {
                                    radius = MIN_ZOOM;
                                    polar_angle = 0.;
                                    azimuthal_angle = 0.;
                                    ui.close_menu();
                                }
                                if ui.button("Front").clicked() {
                                    radius = MAX_ZOOM;
                                    polar_angle = 0.;
                                    azimuthal_angle = 0.;
                                    ui.close_menu();
                                }
                                if ui.button("Left").clicked() {
                                    radius = MAX_ZOOM;
                                    polar_angle = 0.;
                                    azimuthal_angle = FRAC_PI_2;
                                    ui.close_menu();
                                }
                                if ui.button("Top").clicked() {
                                    radius = MAX_ZOOM;
                                    polar_angle = MAX_POLAR;
                                    azimuthal_angle = 0.;
                                    ui.close_menu();
                                }
                            });
                            egui::ComboBox::new("Color Scheme Box", "")
                                .width(140.)
                                .selected_text(format!("{}", color_scheme))
                                .show_ui(ui, |ui| {
                                    //ui.selectable_value(&mut color_scheme, ColorScheme::Charge, "Charge");
                                    ui.selectable_value(&mut color_scheme, ColorScheme::Phase, ColorScheme::Phase.to_string());
                                    //ui.selectable_value(&mut color_scheme, ColorScheme::Distance, ColorScheme::Distance.to_string());
                                    ui.selectable_value(&mut color_scheme, ColorScheme::Probability, ColorScheme::Probability.to_string());
                                });
                            ui.checkbox(&mut show_bohr_radius, "Show Bohr orbit");
                            ui.checkbox(&mut show_crosshair, "Show crosshair");
                            ui.checkbox(&mut show_debug, "Show debug");
                        });
                    });
                    ui.add(egui::Separator::default().vertical());
                    // Right
                    let mut orbital_popup_layout = egui::text::LayoutJob::default();
                    wf.layout_orbital_name(&mut orbital_popup_layout);
                    let orbital_popup_response = ui.add_sized(egui::Vec2::new(45., 30.), egui::Button::new(orbital_popup_layout));
                    let orbital_popup_id = ui.make_persistent_id("Orbital Popup Button");
                    if orbital_popup_response.clicked() {
                        ui.memory().toggle_popup(orbital_popup_id);
                    }
                    egui::popup_below_widget(ui, orbital_popup_id, &orbital_popup_response, |ui| {
                        ui.set_min_width(100.);
                        let mut n = wf.get_n();
                        ui.add(egui::Slider::new(&mut n, 1..=5).integer().step_by(1.).text("Principal (n)"));
                        let mut l = wf.get_l().clamp(0, n-1);
                        if n > 1 {
                            ui.add(egui::Slider::new(&mut l, 0..=n-1).integer().step_by(1.).text("Azimuthal (ℓ)"));
                        }
                        let mut m = wf.get_m().clamp(-(l as i32), l as i32);
                        if l > 0 {
                            ui.add(egui::Slider::new(&mut m, -(l as i32)..=l as i32).integer().step_by(1.).text("Magnetic (m)"));
                        }
                        if n != wf.get_n() || l != wf.get_l() || m != wf.get_m() {
                            electrons.clear();
                            wf = Wavefunction::new(n, l, m);
                        }
                    });
                    ui.vertical(|ui| {
                        ui.label("Orbital");
                        ui.horizontal(|ui| {
                            // FIXME: Figure out why keyboard input doesn't work
                            ui.add(egui::Slider::new(&mut num_particles, 1..=5_000).integer().show_value(true).text("x"));
                            if ui.button("Measure").clicked() {
                                spawn_electrons(&wf, num_particles, &mut electrons, &mut max_prob, &mut max_rad);
                            }
                        });
                    })
                });
            });
            /*egui::Window::new("Config").default_pos(Pos2::new(100., 100.)).show(ctx, |ui| {
                ui.label("Render");
            });*/
            egui_using_pointer = ctx.is_using_pointer() || ctx.is_pointer_over_area();
        });
        egui_macroquad::draw();

        next_frame().await
    }
}

fn spawn_electrons(wf: &Wavefunction, num: i32, electrons: &mut Vec<Electron>, max_prob: &mut f32, max_rad: &mut f32) {
    electrons.clear();
    *max_prob = 0.;
    *max_rad = 0.;
    for _ in 0..num {
        let e = wf.measure();
        if e.total_prob > *max_prob {
            *max_prob = e.total_prob;
        }
        if e.r > *max_rad {
            *max_rad = e.r;
        }
        electrons.push(e);
    }
}

// Number of sides to draw for particles, or 0 if particle should be culled
fn particle_lod(camera_pos: Vec3, particle_pos: Vec3, forward: Vec3, cull_fov: f32) -> u8 {
    if (particle_pos - camera_pos).angle_between(forward).abs() > cull_fov.to_radians() {
        return 0
    }
    let dist_squared = particle_pos.distance_squared(camera_pos);
    if dist_squared < 100. {
        20 - (0.14 * dist_squared).round() as u8
    } else {
        6
    }
}

type Vertex = ([f32; 3], [f32; 2], [f32; 4]);
fn new_vertex(x: f32, y: f32, z: f32, u: f32, v: f32, color: Color) -> Vertex {
    /* up (0,1,0) right (1,0,0)
    x*(1,0,0) + y*(0,1,0)
    (x,0,0) + (0,y,0)
     */
    ([x, y, z], [u, v], [color.r, color.g, color.b, color.a])
}

fn draw_poly_billboard(quad_gl: &mut QuadGl, up: Vec3, right: Vec3, pos: Vec3, sides: u8, radius: f32, rotation: f32, color: Color) {
    if sides == 0 {
        return;
    }
    let mut vertices = Vec::<Vertex>::with_capacity(sides as usize + 2);
    let mut indices = Vec::<u16>::with_capacity(sides as usize * 3);

    let face = vec3(1.,1.,1.).normalize();//up.add(right).normalize();
    let rot = rotation.to_radians();
    vertices.push(new_vertex(pos.x, pos.y, pos.z, 0., 0., color));
    for i in 0..sides + 1 {
        let rx = (i as f32 / sides as f32 * std::f32::consts::PI * 2. + rot).cos();
        let ry = (i as f32 / sides as f32 * std::f32::consts::PI * 2. + rot).sin();
        let v_pos = pos + up.mul(ry*radius) + right.mul(rx*radius);
        let vertex = new_vertex(v_pos.x, v_pos.y, v_pos.z, rx, ry, color);

        vertices.push(vertex);

        if i != sides {
            indices.extend_from_slice(&[0, i as u16 + 1, i as u16 + 2]);
        }
    }

    quad_gl.texture(None);
    quad_gl.draw_mode(DrawMode::Triangles);
    quad_gl.geometry(&vertices, &indices);
}