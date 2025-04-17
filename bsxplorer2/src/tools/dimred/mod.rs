pub mod dist;
pub mod mle;
pub mod mom;

use num::NumCast;

#[allow(unused)]
struct MethCountsDataOwned {
    positions: Vec<u32>,
    count_m: Vec<u16>,
    count_total: Vec<u16>,
    density: Vec<f64>,
    density_cumsum: Vec<f64>,
}

fn to_num<F, T>(x: F) -> T
where
    T: NumCast,
    F: num::ToPrimitive,
{
    T::from(x).unwrap()
}
