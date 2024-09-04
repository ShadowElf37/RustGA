use std::fmt::Display;
use std::ops::{self, Range};

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



pub struct GeometricAlgebra {
    dimension: usize
}
impl GeometricAlgebra {
    const MAX_DIMENSION: usize = 12; // size of multivector will be 2^MAX_DIMENSION bytes so please be very careful changing this!

    pub fn new(dimension: usize) -> Self {
        assert!(dimension <= Self::MAX_DIMENSION, "GA dimension must be <= {} due to memory constraints", Self::MAX_DIMENSION);
        Self {dimension: dimension}
    }

    pub fn zero(&self) -> Multivector {
        Multivector::new_with_dimension(self.dimension)
    }
    pub fn e(&self, i: usize) -> Multivector {
        assert!(i <= self.dimension, "Cannot make e{} in a {}-dimensional GA", i, self.dimension);
        let mut m = self.zero();
        m.blades[1 << i-1] = 1.0;
        m
    }
}


impl ops::AddAssign for Multivector {
    fn add_assign(&mut self, rhs: Self) {
        for i in 0..self.blades.len() {
            self.blades[i] += rhs.blades[i];
        }
    }
}
// no mul assign because it has to do extra allocations regardless, just use a = a * b

impl ops::Add<Multivector> for Multivector {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        let mut new = self.clone();
        new += rhs;
        new
    }
}
impl ops::Mul<Multivector> for Multivector {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut new = Multivector::new_like(&self);
        //let my_nonzero = ;
        //let their_nonzero = ;
        for my_blade in (0..self.blades.len()).filter(|i| self.blades[*i] != 0.0) {
            for their_blade in (0..rhs.blades.len()).filter(|i| rhs.blades[*i] != 0.0) {
                let scalar1 = self.blades[my_blade];
                let scalar2 = rhs.blades[their_blade];
                if scalar1 != 0.0 && scalar2 != 0.0 {
                    let (elements, sign) = self.multiply_two_normalized_blades(my_blade as Blade, their_blade as Blade);
                    new.blades[elements as usize] += scalar1 * scalar2 * sign;
                }
            }
        }
        new
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
        for my_blade in 0..self.blades.len() {
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

    // pub fn dot(self, other: Multivector) -> Self {
    //     return (self * other).project_out_grade(grade)
    // }

    // returns new blade and the sign
    fn multiply_two_normalized_blades(&self, lhs: Blade, rhs: Blade) -> (Blade, Scalar) {
        let mut permutation_sign = 0;
        if lhs != 0 && rhs != 0 {
            let mut jumps = 0;
            for i in (0..=highest_bit(lhs)).rev() {
                permutation_sign += jumps * get_bit(rhs, i);
                jumps += get_bit(lhs, i);
            }
        }

        (lhs ^ rhs, [1.0, -1.0][(permutation_sign % 2) as usize])
    }
}



fn main() {
    let g = GeometricAlgebra::new(3);
    let a = g.e(1) + g.e(2);
    let b = g.e(2) * g.e(1);
    println!("{:?}", -32.1*g.e(3)*(a.clone()+b.clone()));
    println!("{}", a.clone()*b.clone());
    println!("{} {}", a, b);
}
