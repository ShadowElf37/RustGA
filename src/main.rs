use std::fmt::Display;
use std::f64::consts::*;

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

impl Multivector {
    fn nonzero_blade_indices(&self) -> impl Iterator<Item = Blade> + '_ {
        (0..self.blades.len()).filter(|i| self.blades[*i] != 0.0)
    }
}

pub struct GeometricAlgebra {
    pub dimension: usize,
    pub signature: (u32, u32),
    p_mask: usize,
    q_mask: usize,
    //signs_of_squares: Vec<f64>,
    multiplication_table: Vec<u8>,
    size: usize,
}
impl GeometricAlgebra {
    const MAX_DIMENSION: u32 = 12;
    // size of multivector will be 2^MAX_DIMENSION bytes
    // size of the internal multiplication table will be 4^MAX_DIMENSION bytes
    // please be very careful changing this! your computer might explode!

    pub fn new(dimension: u32, signature: (u32, u32)) -> Self {
        assert!(dimension <= Self::MAX_DIMENSION, "GA dimension must be <= {} due to memory constraints. You can override this, but please be very careful.", Self::MAX_DIMENSION);
        assert!(signature.0 + signature.1 == dimension, "Signature p+q do not add up to the dimension of the algebra");
        let q_mask = usize::MAX << signature.0;
        let p_mask = !q_mask;
        let size = 2usize.pow(dimension);

        let multiplication_table = {
            let mut v = Vec::<u8>::with_capacity(size*size);
            for lhs in 0..size {
                for rhs in 0..size {
                    let mut permutation_sign = (lhs & rhs & q_mask).count_ones() as usize;
                    if lhs != 0 && rhs != 0 {
                        let mut jumps = 0;
                        for i in (0..=highest_bit(lhs)).rev() {
                            permutation_sign += jumps * get_bit(rhs, i);
                            jumps += get_bit(lhs, i);
                        }
                    }
                    v.push((permutation_sign & 1) as u8);
                }
            }
            v
        };

        println!("Geometric algebra loaded with dimension {} and signature {:?}. Multiplication table is {} bytes.", dimension, signature, size*size);
        /*let signs_of_squares = {
            let mut v = Vec::with_capacity(size);
            for blade in 0..size {
                let p = (blade & p_mask).count_ones();
                let q = (blade & q_mask).count_ones();
                v.push([1.0, 1.0, -1.0, -1.0][((p as i32 - q as i32).rem_euclid(4)) as usize]);
            }
            v
        };*/

        Self {
            dimension: dimension as usize,
            signature,
            p_mask,
            q_mask,
            multiplication_table,
            size
        }
    }

    pub fn exp(&self, m: &Multivector) -> Multivector {
        let mut result = self.scalar(1.0);
        for blade in m.nonzero_blade_indices() {
            let val = m.blades[blade];
            match self.sign_of_square(blade) {
                1.0 => {
                    let r1 = self.scale(&result, val.cosh());
                    let r2 = self.mul_blade_R(&result, (val.sinh(), blade));
                    result = self.add(&r1, &r2);
                }
                -1.0 => {
                    let r1 = self.scale(&result, val.cos());
                    let r2 = self.mul_blade_R(&result, (val.sin(), blade));
                    result = self.add(&r1, &r2);
                } 
                _ => panic!()
            }
        }
        result
    }

    pub fn log(&self, m: &Multivector) -> Multivector {
        let mut result = m.clone();

        let mut several_blades = false;
        let mut the_only_blade: Blade = 0;
        for blade in 1..self.size {
            if m.blades[blade] != 0.0 {
                if the_only_blade == 0 {
                    the_only_blade = blade;
                }
                else {
                    several_blades = true;
                    break;
                }
            }
        }


        if !several_blades{
            if the_only_blade == 0 {
                result.blades[0] = result.blades[0].ln();
                return result;
            } else {
                match self.sign_of_square(the_only_blade) {
                    1.0 => {
                        result.blades[0] = self.quadratic_form(m).sqrt().ln();
                        result.blades[the_only_blade] = Scalar::atanh(m.blades[the_only_blade]/m.blades[0]);
                    }
                    -1.0 => {
                        result.blades[0] = self.quadratic_form(m).sqrt().ln();
                        result.blades[the_only_blade] = Scalar::atan(m.blades[the_only_blade]/m.blades[0]);
                    }
                    _ => panic!()
                }
                return result;
            }
            
        }
        else { // compute iteratively using Newton's method
            panic!("cannot compute logarithms of multivectors with multiple non-scalar blades");
            let mag = self.sum_of_squares(&result).sqrt();
            result = self.scale(&result, 1.0/mag);
            result = self.add_scalar(&result, -1.0);

            let initial = result.clone();
            let mut product = result.clone();

            for n in 2..30 {
                product = self.mul(&product, &initial);
                result = self.add(&result, &self.scale(&product, [-1.0, 1.0][n&1] / (n as f64)));

                // let exp = self.exp(&result);
                // result = self.add(&result, &self.scale(
                //     &self.mul(&self.sub(m, &exp), &self.inverse(&self.add(m, &exp))), 2.0
                // ));
            }
            return self.add_scalar(&result, mag.ln());
        } 
    }

    // MULTIPLICATION HELPERS
        fn sign_of_square(&self, blade: Blade) -> f64 {
            self.get_sign_of_multiplication_of_two_blades(blade, blade)
        }
        fn mul_blade_R(&self, lhs: &Multivector, rhs: (Scalar, Blade)) -> Multivector {
            let mut new = self.zero();
            for blade in lhs.nonzero_blade_indices() {
                let elements = blade ^ rhs.1;
                new.blades[elements as usize] += lhs.blades[blade] * rhs.0 * self.get_sign_of_multiplication_of_two_blades(blade, rhs.1);
            }
            new
        }
        fn mul_blade_L(&self, lhs: (Scalar, Blade), rhs: &Multivector) -> Multivector {
            let mut new = self.zero();
            for blade in rhs.nonzero_blade_indices() {
                let elements = blade ^ lhs.1;
                new.blades[elements as usize] += rhs.blades[blade] * lhs.0 * self.get_sign_of_multiplication_of_two_blades(lhs.1, blade);
            }
            new
        }
        fn get_sign_of_multiplication_of_two_blades_manually(&self, lhs: Blade, rhs: Blade) -> Scalar {
            // for small dimension, this could be a table. consider - at dim=12, the table is only 2^24 large (16 MB).
            // fall back on this if the code is modified to support larger dimensions
            let mut permutation_sign = (lhs & rhs & self.q_mask).count_ones() as usize;
            if lhs != 0 && rhs != 0 {
                let mut jumps = 0;
                for i in (0..=highest_bit(lhs)).rev() {
                    permutation_sign += jumps * get_bit(rhs, i);
                    jumps += get_bit(lhs, i);
                }
            }
            [1.0, -1.0][(permutation_sign % 2) as usize]
        }
        fn get_sign_of_multiplication_of_two_blades(&self, lhs: Blade, rhs: Blade) -> Scalar {
            [1.0, -1.0][self.fetch_multiplication_table(lhs, rhs) as usize]
        }
        fn fetch_multiplication_table(&self, lhs: Blade, rhs: Blade) -> u8 {
            self.multiplication_table[lhs*self.size + rhs]
        }

    pub fn add(&self, lhs: &Multivector, rhs: &Multivector) -> Multivector {
        let mut new = lhs.clone();
        for blade in 0..self.size {
            new.blades[blade] += rhs.blades[blade]
        }
        new
    }
    pub fn add_scalar(&self, lhs: &Multivector, rhs: Scalar) -> Multivector {
        let mut new = lhs.clone();
        new.blades[0] += rhs;
        new
    }
    pub fn sub(&self, lhs: &Multivector, rhs: &Multivector) -> Multivector {
        let mut new = lhs.clone();
        for blade in 0..self.size {
            new.blades[blade] -= rhs.blades[blade]
        }
        new
    }
    pub fn neg(&self, m: &Multivector) -> Multivector {
        let mut new = m.clone();
        for blade in 0..self.size {
            new.blades[blade] = -new.blades[blade]
        }
        new
    }
    pub fn scale(&self, lhs: &Multivector, rhs: Scalar) -> Multivector {
        let mut new = lhs.clone();
        for blade in 0..self.size {
            new.blades[blade] *= rhs;
        }
        new
    }
    pub fn mul(&self, lhs: &Multivector, rhs: &Multivector) -> Multivector {
        let mut new = self.zero();
        for my_blade in lhs.nonzero_blade_indices() {
            for their_blade in rhs.nonzero_blade_indices() {
                let elements = my_blade ^ their_blade;
                new.blades[elements as usize] += lhs.blades[my_blade] * rhs.blades[their_blade] * self.get_sign_of_multiplication_of_two_blades(my_blade, their_blade);
            }
        }
        new
    }
    pub fn mul_consecutive(&self, multivectors: Vec<&Multivector>) -> Multivector {
        let mut result = self.scalar(1.0);
        for m in multivectors {
            result = self.mul(&result, m);
        }
        result
    }

    pub fn commutator_product(&self, lhs: &Multivector, rhs: &Multivector) -> Multivector {
        self.scale(&self.sub(&self.mul(lhs, rhs), &self.mul(rhs, lhs)), 0.5)
    }
    pub fn inner_product(&self, lhs: &Multivector, rhs: &Multivector) -> Multivector {
        let mut new = self.zero();
        for my_blade in lhs.nonzero_blade_indices() {
            let my_grade = my_blade.count_ones() as i32;
            for their_blade in rhs.nonzero_blade_indices() {
                let elements = my_blade ^ their_blade;
                if elements.count_ones() as i32 == (my_grade - their_blade.count_ones() as i32).abs() {
                    new.blades[elements as usize] += lhs.blades[my_blade] * rhs.blades[their_blade] * self.get_sign_of_multiplication_of_two_blades(my_blade, their_blade);
                }
            }
        }
        new
    }
    pub fn outer_product(&self, lhs: &Multivector,  rhs: &Multivector) -> Multivector {
        let mut new = self.zero();
        for my_blade in lhs.nonzero_blade_indices() {
            for their_blade in rhs.nonzero_blade_indices() {
                if their_blade & my_blade == 0 { // if they have none in common, since this is the antisymmetric product
                    new.blades[(their_blade ^ my_blade) as usize] += lhs.blades[my_blade] * rhs.blades[their_blade] * self.get_sign_of_multiplication_of_two_blades(my_blade, their_blade);
                }
            }
        }
        new
    }

    pub fn project_out_grade(&self, m: &Multivector, grade: u32) -> Multivector {
        let mut new = self.zero();
        for blade in 0..self.size {
            if blade.count_ones() == grade {
                new.blades[blade] = m.blades[blade]
            }
        }
        new
    }

    // returns None if it's not homogeneous, or it returns the only blade it has
    pub fn is_homogeneous(&self, m: &Multivector) -> Option<Blade> {
        let mut grade: Option<u32> = None;
        let mut only_blade: Option<Blade> = None;
        for blade in m.nonzero_blade_indices() {
            if grade.is_none() {
                grade = Some(blade.count_ones());
                only_blade = Some(blade);
                continue;
            }
            if blade.count_ones() != grade.unwrap() {
                return None;
            }
        }
        return Some(only_blade.unwrap_or(0));
    }
    // counts e.g.  1.0 + e1e2  as homogeneous
    pub fn is_homogeneous_except_scalar_part(&self, m: &Multivector) -> Option<Blade> {
        let mut grade: Option<u32> = None;
        let mut only_blade: Option<Blade> = None;
        for blade in m.nonzero_blade_indices() {
            if blade == 0 {
                continue;
            }
            if grade.is_none() {
                grade = Some(blade.count_ones());
                only_blade = Some(blade);
                continue;
            }
            if blade.count_ones() != grade.unwrap() {
                return None;
            }
        }
        return Some(only_blade.unwrap_or(0));
    }
    pub fn is_homogeneous_of_grade(&self, m: &Multivector, grade: u32) -> bool {
        for blade in m.nonzero_blade_indices() {
            if blade.count_ones() != grade {
                return false;
            }
        }
        return true;
    }
    pub fn is_in_center_of_algebra(&self, m: &Multivector) -> bool {
        for blade1 in m.nonzero_blade_indices() {
            for blade2 in 0..self.size {
                //println!("{}, {}, {}, {}", self.fetch_multiplication_table(blade1, blade2), self.fetch_multiplication_table(blade2, blade1), blade1, blade2);
                if self.fetch_multiplication_table(blade1, blade2) != self.fetch_multiplication_table(blade2, blade1) {
                    return false;
                }
            }
        }
        return true;
    }

    pub fn mul_and_project_out_scalar(&self, lhs: &Multivector, rhs: &Multivector) -> Scalar {
        let mut result: Scalar = 0.0;
        for blade in 0..self.size {
            result += lhs.blades[blade] * rhs.blades[blade] * self.sign_of_square(blade)
        }
        result
    }

    // generalization of complex conjugate
    pub fn reverse(&self, m: &Multivector) -> Multivector {
        let mut new = m.clone();
        for blade in m.nonzero_blade_indices() {
            new.blades[blade] *= [1.0, -1.0][((blade.count_ones()/2) % 2) as usize]
        }
        new
    }

    pub fn rotate_in(&self, m: &Multivector, plane: &Multivector, theta: Scalar) -> Multivector {
        let argument = self.scale(plane, -theta / 2.0);
        let r = &self.exp(&argument);
        self.conjugate(m, r)
    }

    // for multivector m and conjugator R, compute RmR*
    pub fn conjugate(&self, m: &Multivector, conjugator: &Multivector) -> Multivector {
        self.mul(conjugator, &self.mul(m, &self.reverse(conjugator)))
    }

    // hodge star
    pub fn dual(&self, m: &Multivector) -> Multivector {
        self.mul_blade_L((1.0, self.size-1), m)
    }

    pub fn quadratic_form(&self, m: &Multivector) -> Scalar {
        let mut sum: Scalar = 0.0;
        for blade in 0..self.size {
            sum += m.blades[blade] * m.blades[blade] * [1.0, -1.0][((blade & self.q_mask).count_ones() & 1) as usize]
        }
        sum
        //self.mul_and_project_out_scalar(m, &self.reversed(m))
    }
    pub fn sum_of_squares(&self, m: &Multivector) -> Scalar {
        let mut sum: Scalar = 0.0;
        for blade in 0..self.size {
            sum += m.blades[blade] * m.blades[blade]
        }
        sum
    }

    pub fn inverse(&self, m: &Multivector) -> Multivector {
        self.scale(m, 1.0/self.quadratic_form(m))
    }


    // HELPER INITIALIZERS
        fn new_mv_uninit(&self) -> Multivector {
            Multivector{blades: Vec::<Scalar>::with_capacity(self.size)}
        }
        fn new_mv_with_value_at_index(&self, index: Blade, value: Scalar) -> Multivector {
            let mut new = self.zero();
            new.blades[index] = value;
            new
        }
        fn new_unit_mv_from_blade(&self, blade: Blade) -> Multivector {
            let mut new = self.zero();
            new.blades[blade] = 1.0;
            new
        }

    // MULTIVECTOR GENERATORS
        pub fn zero(&self) -> Multivector {
            let mut new = self.new_mv_uninit();
            for _ in 0..self.size {
                new.blades.push(0f64)
            }
            new
        }
        pub fn scalar(&self, x: Scalar) -> Multivector {
            self.new_mv_with_value_at_index(0, x)
        }
        pub fn e(&self, i: usize) -> Multivector {
            assert!(i <= self.dimension, "Cannot make e{} in a {}-dimensional GA", i, self.dimension);
            self.new_mv_with_value_at_index(1 << i-1, 1.0)
        }
        pub fn blade(&self, e: &[usize]) -> Multivector {
            let mut new = self.zero();
            let mut blade = 0;
            for i in e {
                assert!(*i <= self.dimension, "Cannot make e{} in a {}-dimensional GA", i, self.dimension);
                blade |= 1 << *i-1;
            }
            new.blades[blade] = 1.0;
            new
        }
        // complex number
        pub fn scalar_pseudoscalar(&self, a: Scalar, b: Scalar) -> Multivector {
            assert!(self.dimension >= 2, "Cannot make complex number in a {}-dimensional GA", self.dimension);
            let mut new = self.pseudoscalar(b);
            new.blades[0] = a;
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
        pub fn vec4(&self, t: Scalar, x: Scalar, y: Scalar, z: Scalar) -> Multivector {
            assert!(self.dimension >= 4, "Cannot make vec4 in a {}-dimensional GA", self.dimension);
            let mut new = self.zero();
            new.blades[0b0001] = t;
            new.blades[0b0010] = x;
            new.blades[0b0100] = y;
            new.blades[0b1000] = z;
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
        pub fn bivec4(&self, xt: Scalar, yt: Scalar, zt: Scalar, xy: Scalar, yz: Scalar, zx: Scalar) -> Multivector {
            assert!(self.dimension >= 4, "Cannot make 6-element bivector in a {}-dimensional GA", self.dimension);
            let mut new = self.zero();
            new.blades[0b0011] = xt;
            new.blades[0b0101] = yt;
            new.blades[0b1001] = zt;
            new.blades[0b0110] = xy;
            new.blades[0b1100] = yz;
            new.blades[0b1010] = zx;
            new
        }
        pub fn trivec4(&self, xyz: Scalar, tyz: Scalar, tzx: Scalar, txy: Scalar) -> Multivector {
            assert!(self.dimension >= 4, "Cannot make 4-element trivector in a {}-dimensional GA", self.dimension);
            let mut new = self.zero();
            new.blades[0b1110] = xyz;
            new.blades[0b1101] = tyz;
            new.blades[0b1011] = tzx;
            new.blades[0b0111] = txy;
            new
        }
        pub fn pseudoscalar(&self, x: Scalar) -> Multivector {
            self.new_mv_with_value_at_index(self.size-1, x)
        }
}

fn main() {
    let g = GeometricAlgebra::new(4, (1, 3));
    let a = &g.vec2(1.0, 1.0);
    let c = &g.vec2(1.0, 4.0);
    println!("{} {} {}", g.dual(&g.e(1)), g.dual(&g.e(2)), g.dual(&g.e(3)));
    println!("{}", g.outer_product(c, a));
    println!("{}", g.rotate_in(a, &g.blade(&[1, 2]), PI/4.0));
    println!("{}", g.exp(&g.scalar(1.0)));
    println!("{}", g.add(&g.scalar(1.0), &g.pseudoscalar(2.0)));
    println!("{}", g.quadratic_form(&g.vec4(2.0, 0.0, 1.0, 1.0)));
    println!("{}", g.quadratic_form(&g.vec4(2.0, 0.0, 1.0, 1.0)));

    let g = GeometricAlgebra::new(2, (2, 0));
    let angle = g.scalar_pseudoscalar(2.0, PI/4.0);
    let angle = g.add(&angle, &g.e(1));
    let n = g.exp(&angle);
    println!("{}", angle);
    println!("{}", n);
    println!("{:?}", g.multiplication_table);
    

    let g = GeometricAlgebra::new(6, (6, 0));
}
