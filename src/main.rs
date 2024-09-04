use std::fmt::Display;
use std::ops::{self, Range};
use std::f64::consts::*;

static FRACTIONS: [f64; 64] = [1.0, 1.0, 0.5, 0.3333333333333333, 0.25, 0.2, 0.16666666666666666, 0.14285714285714285, 0.125, 0.1111111111111111, 0.1, 0.09090909090909091, 0.08333333333333333, 0.07692307692307693, 0.07142857142857142, 0.06666666666666667, 0.0625, 0.058823529411764705, 0.05555555555555555, 0.05263157894736842, 0.05, 0.047619047619047616, 0.045454545454545456, 0.043478260869565216, 0.041666666666666664, 0.04, 0.038461538461538464, 0.037037037037037035, 0.03571428571428571, 0.034482758620689655, 0.03333333333333333, 0.03225806451612903, 0.03125, 0.030303030303030304, 0.029411764705882353, 0.02857142857142857, 0.027777777777777776, 0.02702702702702703, 0.02631578947368421, 0.02564102564102564, 0.025, 0.024390243902439025, 0.023809523809523808, 0.023255813953488372, 0.022727272727272728, 0.022222222222222223, 0.021739130434782608, 0.02127659574468085, 0.020833333333333332, 0.02040816326530612, 0.02, 0.0196078431372549, 0.019230769230769232, 0.018867924528301886, 0.018518518518518517, 0.01818181818181818, 0.017857142857142856, 0.017543859649122806, 0.017241379310344827, 0.01694915254237288, 0.016666666666666666, 0.01639344262295082, 0.016129032258064516, 0.015873015873015872];

type Blade = usize;
type Scalar = f64;

fn get_bit(a: Blade, bit_index: usize) -> usize {
    ((1 << bit_index) & a) >> bit_index
}
fn highest_bit(a: Blade) -> usize {
    a.ilog2() as usize
}


#[derive(Debug, Clone)]
pub struct Multivector{
    blades: Vec<Scalar>
}

impl Display for Multivector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut printed_any_vectors = false;
        write!(f, "Multivector(")?;
        for i in 0..self.blades.len() {

            let scalar = self.blades[i];

            if scalar != 0.0 {
                if printed_any_vectors {
                    if scalar < 0.0 {
                        write!(f, " - ")?;
                    } else {
                        write!(f, " + ")?;
                    }
                } else {
                    if scalar < 0.0 {
                        write!(f, "-")?;
                    }
                }
                printed_any_vectors = true;

                if i != 0 {
                    if scalar.abs() != 1.0 {
                        write!(f, "{}*", scalar.abs())?;
                    }

                    for j in 0..=highest_bit(i as Blade) {
                        match get_bit(i as Blade, j) {
                            1 => write!(f, "e{}", j+1)?,
                            _ => continue,
                        };
                    }
                } else {
                    write!(f, "{}", self.blades[i].abs())?;
                }
                
            }            
        }
        if !printed_any_vectors {
            write!(f, "{}", self.blades[0])?;
        }
        write!(f, ")")?;
        Ok(())
    }
}

impl ops::Add<Multivector> for Multivector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut new = self.clone();
        new += rhs;
        new
    }
}
impl ops::AddAssign for Multivector {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..self.blades.len() {
            self.blades[i] += rhs.blades[i];
        }
    }
}
impl ops::Mul<Multivector> for Multivector {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut new = Multivector::new_like(&self);
        for my_blade in self.nonzero_blade_indices() {
            for their_blade in rhs.nonzero_blade_indices() {
                let elements = my_blade ^ their_blade;
                new.blades[elements as usize] += self.blades[my_blade] * rhs.blades[their_blade] * Self::get_sign_of_multiplication_of_two_blades(my_blade, their_blade);
            }
        }
        new
    }
}
impl ops::MulAssign<Multivector> for Multivector {
    fn mul_assign(&mut self, rhs: Self) {
        self.blades = (self.clone() * rhs).blades;
    }
}

impl ops::Add<Scalar> for Multivector {
    type Output = Self;

    fn add(self, rhs: Scalar) -> Self::Output {
        let mut new = self.clone();
        new += rhs;
        new
    }
}
impl ops::AddAssign<Scalar> for Multivector {
    fn add_assign(&mut self, rhs: Scalar) {
        self.blades[0] += rhs
    }
}
impl ops::Mul<Scalar> for Multivector {
    type Output = Self;

    fn mul(self, rhs: Scalar) -> Self::Output {
        let mut new = self.clone();
        new *= rhs;
        new
    }
}
impl ops::MulAssign<Scalar> for Multivector {
    fn mul_assign(&mut self, rhs: Scalar) {
        for my_blade in self.blade_indices() {
            self.blades[my_blade] *= rhs
        }
    }
}
impl ops::Mul<Multivector> for Scalar {
    type Output = Multivector;

    fn mul(self, rhs: Multivector) -> Self::Output {
        let mut new = rhs.clone();
        new *= self;
        new
    }
}
impl ops::Add<Multivector> for Scalar {
    type Output = Multivector;

    fn add(self, rhs: Multivector) -> Self::Output {
        let mut new = rhs.clone();
        new += self;
        new
    }
}

impl Multivector {
    fn new_with_size(size: usize) -> Self {
        Self {
            blades: (0..size).map(|_| 0f64).collect::<Vec<Scalar>>()
        }
    }
    fn new_with_dimension(dimension: usize) -> Self {
        Self::new_with_size(2usize.pow(dimension as u32))
    }
    fn new_like(other: &Multivector) -> Self {
        Self::new_with_size(other.blades.len())
    }

    fn blade_indices(&self) -> Range<Blade> {
        0..self.blades.len()
    }
    fn nonzero_blade_indices(&self) -> impl Iterator<Item = Blade> + '_ {
        (0..self.blades.len()).filter(|i| self.blades[*i] != 0.0)
    }
    
    pub fn dimension(&self) -> usize {
        self.blades.len().ilog2() as usize
    }
    pub fn project_out_grade(&self, grade: u32) -> Self {
        let mut new = Multivector::new_like(&self);
        for blade in self.blade_indices() {
            if blade.count_ones() == grade {
                new.blades[blade] = self.blades[blade]
            }
        }
        new
    }

    pub fn is_homogeneous(&self) -> bool {
        let mut grade: Option<u32> = None;
        for blade in self.nonzero_blade_indices() {
            if grade.is_none() {
                grade = Some(blade.count_ones());
                continue;
            }
            if blade.count_ones() != grade.unwrap() {
                return false;
            }
        }
        return true;
    }
    pub fn is_homogeneous_of_grade(&self, grade: u32) -> bool {
        for blade in self.nonzero_blade_indices() {
            if blade.count_ones() != grade {
                return false;
            }
        }
        return true;
    }

    fn get_sign_of_multiplication_of_two_blades(lhs: Blade, rhs: Blade) -> Scalar {
        let mut permutation_sign = 0;
        if lhs != 0 && rhs != 0 {
            let mut jumps = 0;
            for i in (0..=highest_bit(lhs)).rev() {
                permutation_sign += jumps * get_bit(rhs, i);
                jumps += get_bit(lhs, i);
            }
        }
        [1.0, -1.0][(permutation_sign % 2) as usize]
    }

    pub fn multiply_and_project_out_scalar(&self, rhs: &Multivector) -> Scalar {
        let mut result: Scalar = 0.0;
        for blade in self.blade_indices() {
            result += self.blades[blade] * rhs.blades[blade]
        }
        result
    }

    pub fn commutator_product(&self, rhs: &Multivector) -> Multivector {
        0.5 * (self.clone() * rhs.clone() + rhs.clone() * self.clone())
    }

    pub fn inner_product(&self, rhs: &Multivector) -> Multivector {
        let mut new = Multivector::new_like(&self);
        for my_blade in self.nonzero_blade_indices() {
            for their_blade in rhs.nonzero_blade_indices() {
                let elements = my_blade ^ their_blade;
                if elements.count_ones() as i32 == (my_blade.count_ones() as i32 - their_blade.count_ones() as i32).abs() {
                    new.blades[elements as usize] += self.blades[my_blade] * rhs.blades[their_blade] * Self::get_sign_of_multiplication_of_two_blades(my_blade, their_blade);
                }
            }
        }
        new
    }
    pub fn outer_product(&self, rhs: &Multivector) -> Multivector {
        let mut new = Multivector::new_like(&self);
        for my_blade in self.nonzero_blade_indices() {
            for their_blade in rhs.nonzero_blade_indices() {
                if their_blade & my_blade == 0 { // if they have none in common, since this is the antisymmetric product
                    new.blades[(their_blade ^ my_blade) as usize] += self.blades[my_blade] * rhs.blades[their_blade] * Self::get_sign_of_multiplication_of_two_blades(my_blade, their_blade);
                }
            }
        }
        new
    }

    pub fn mag2(&self) -> Scalar {
        self.blades.iter().map(|x| x*x).sum()
    }

    pub fn rotate_in(self, plane: &Multivector, theta: Scalar) -> Multivector {
        let argument = plane.clone() * (theta / 2.0);
        let r = (argument.clone() * -1.0).exp();
        let rdagger = argument.clone().exp();
        r * self * rdagger
    }

    pub fn reverse(&mut self) -> &mut Self {
        for blade in 0..self.blades.len() {
            if self.blades[blade] == 0.0 {continue}
            let r = blade.count_ones();
            let permutation_sign = r*(r.wrapping_sub(1))/2;
            self.blades[blade] *= [1.0, -1.0][(permutation_sign % 2) as usize]
        }
        self
    }

    pub fn exp(&self) -> Multivector {
        let mut product = Self::new_like(&self);
        product.blades[0] = 1.0;
        let mut sum = Self::new_like(&self);
        let taylor_series_anchor_point = self.blades[0].round();
        let e_approx = E.powi(taylor_series_anchor_point as i32);
        sum.blades[0] = 1.0;
        for i in 1..23 {
            product *= (self.clone() + -self.blades[0].round()) * FRACTIONS[i];
            sum += product.clone();
        }
        sum * e_approx
    }
}

pub struct GeometricAlgebra {
    dimension: usize
}
impl GeometricAlgebra {
    const MAX_DIMENSION: usize = 12; // size of multivector will be 2^MAX_DIMENSION bytes so please be very careful changing this!

    pub fn new(dimension: usize) -> Self {
        assert!(dimension <= Self::MAX_DIMENSION, "GA dimension must be <= {} due to memory constraints", Self::MAX_DIMENSION);
        Self {dimension}
    }


    pub fn scalar(&self, x: Scalar) -> Multivector {
        let mut new = self.zero();
        new.blades[0] = x;
        new
    }
    pub fn zero(&self) -> Multivector {
        Multivector::new_with_dimension(self.dimension)
    }
    pub fn e(&self, i: usize) -> Multivector {
        assert!(i <= self.dimension, "Cannot make e{} in a {}-dimensional GA", i, self.dimension);
        let mut new = self.zero();
        new.blades[1 << i-1] = 1.0;
        new
    }
    pub fn pseudoscalar(&self) -> Multivector {
        let mut new = self.zero();
        let last_i = new.blades.len()-1;
        new.blades[last_i] = 1.0;
        new
    }

    pub fn vec2(&self, x: Scalar, y: Scalar) -> Multivector {
        assert!(self.dimension >= 2, "Cannot make vec2 in a {}-dimensional GA", self.dimension);
        let mut new = self.zero();
        new.blades[0b01] = x;
        new.blades[0b10] = y;
        new
    }
    pub fn vec3(&self, x: Scalar, y: Scalar, z: Scalar) -> Multivector {
        assert!(self.dimension >= 3, "Cannot make vec3 in a {}-dimensional GA", self.dimension);
        let mut new = self.zero();
        new.blades[0b001] = x;
        new.blades[0b010] = y;
        new.blades[0b100] = z;
        new
    }
    pub fn vec4(&self, x: Scalar, y: Scalar, z: Scalar, w: Scalar) -> Multivector {
        assert!(self.dimension >= 4, "Cannot make vec4 in a {}-dimensional GA", self.dimension);
        let mut new = self.zero();
        new.blades[0b0001] = x;
        new.blades[0b0010] = y;
        new.blades[0b0100] = z;
        new.blades[0b1000] = w;
        new
    }

    pub fn bivec3(&self, xy: Scalar, yz: Scalar, zx: Scalar) -> Multivector {
        assert!(self.dimension >= 3, "Cannot make 3-element bivector in a {}-dimensional GA", self.dimension);
        let mut new = self.zero();
        new.blades[0b011] = xy;
        new.blades[0b110] = yz;
        new.blades[0b101] = zx;
        new
    }
}

fn main() {
    let g = GeometricAlgebra::new(3);
    let a = g.vec2(1.0, 1.0);
    let c = g.vec2(1.0, 4.0);
    println!("{}", a.inner_product(&c));
    println!("{}", c.outer_product(&a));
    println!("{}", a.rotate_in(&g.pseudoscalar(), PI/4.0));
}
