use macroquad::rand::{rand, gen_range};
use std::cell::{RefMut, RefCell};
use crate::constants::MAX_RADIUS;

/*fn distributed_sample<PDF>(pdf: PDF, pdf_max: f32) -> f32 where PDF: Fn(f32) -> f32 {
    //let mut prob_dist = Uniform::new(0., 1.);
    //let mut sample_dist = Uniform::new(0., 5.*A_0);
    loop {
        let uniform_prob = gen_range(0., 1.);//rng.sample(prob_dist);
        let sample = gen_range(0., 5.);//rng.sample(sample_dist);
        let sample_density = 1.;
        if uniform_prob*pdf_max*sample_density < pdf(sample) {
            return sample;
        }
    }
}*/

pub fn rejection_sample<PDF>(pdf: PDF, pdf_max: f32) -> f32 where PDF: Fn(f32) -> f32 {
    loop {
        let sample_pos = gen_range(0., MAX_RADIUS);
        let sample_prob = pdf(sample_pos);
        let comp_prob = gen_range(0., pdf_max);
        if comp_prob < sample_prob {
            return sample_pos
        }
    }
}

/*const METROPOLIS_WARMUP_ITERS: usize = 10;

pub struct MetropolisSamplerState {
    iters: usize,
    last_pos: f32
}
impl MetropolisSamplerState {
    pub fn new(first_pos: f32) -> Self {
        Self {
            iters: 0,
            last_pos: first_pos
        }
    }
}

fn metropolis_sample_no_warmup<PDF>(state: &RefCell<MetropolisSamplerState>, pdf: &PDF) -> f32 where PDF: Fn(f32) -> f32 {
    let mut mut_state = state.borrow_mut();
    (*mut_state).iters += 1;
    let last_prob: f32 = pdf(mut_state.last_pos);
    let next_pos: f32 = gen_range(0., 3.);
    let next_prob: f32 = pdf(next_pos);

    if next_prob > last_prob {
        (*mut_state).last_pos = next_pos;
    } else {
        let prob_ratio = next_prob / last_prob;
        if (rand() as f32 / u32::MAX as f32) < prob_ratio {
            (*mut_state).last_pos = next_pos;
        }
    }
    mut_state.last_pos
}

pub fn metropolis_sample<PDF>(state: &RefCell<MetropolisSamplerState>, pdf: PDF) -> f32 where PDF: Fn(f32) -> f32 {
    let iters = state.borrow().iters;
    if iters < METROPOLIS_WARMUP_ITERS {
        for _ in iters..METROPOLIS_WARMUP_ITERS {
            metropolis_sample_no_warmup(state, &pdf);
        }
    }
    metropolis_sample_no_warmup(state, &pdf)
}

const CUMUL_INTEGRAL_N: u32 = 100;
const CUMUL_SAMPLE_EPSILON: f32 = 0.1;
pub fn cumulative_sample<PDF>(pdf: PDF, start_pos: f32) -> f32 where PDF: Fn(f32) -> f32 {
    let mut pos = start_pos;
    let mut min = 0.;
    let mut max = f32::MAX;
    loop {
        let cumul_prob = approx_integral(&pdf, 0., pos, CUMUL_INTEGRAL_N);
        if (rand() as f32 / u32::MAX as f32) < cumul_prob {
            max = pos;
        } else {
            min = pos;
        }
        pos = gen_range(min, max);
        if (max - min) < CUMUL_SAMPLE_EPSILON {
            return pos
        }
    }
}*/