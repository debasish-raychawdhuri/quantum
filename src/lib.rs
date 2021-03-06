mod algebra;
use crate::algebra::pow2;
use algebra::*;
use linearalgebra::field::ComplexField;
use linearalgebra::Matrix;
use num_complex::Complex;
extern crate rand;
use rand::prelude::*;

#[derive(PartialEq, Clone, Debug)]
pub struct QCS {
    total_qbits: u8,
    matrix: Matrix<ComplexField>,
}

impl From<Matrix<ComplexField>> for QCS{

    fn from(matrix: Matrix<ComplexField>) -> Self {
        let mut total_qbits = 0;
        let mut dim = matrix.rows();
        loop {
            if dim <= 1 {
                break;
            }
            total_qbits += 1;
            dim >>=1;
        }
        QCS{
            matrix,
            total_qbits
        }
    }
}

impl QCS {
    pub fn new(total_qbits: u8) -> Self {
        let complex_field = ComplexField;
        let matrix = create_swap_matrix(total_qbits, 0, 0, complex_field);
        QCS {
            total_qbits,
            matrix,
        }
    }


    pub fn gate(&self, gate: &QCS, wires: &[u8]) -> QCS {
        let matrix = apply_matrix_to(self.total_qbits, &gate.matrix, wires);
        QCS {
            total_qbits: self.total_qbits,
            matrix: matrix.mul(&self.matrix).unwrap(),
        }
    }

    pub fn controled_gate(&self, gate: &QCS, wires: &[u8], control_wires:&[u8]) -> QCS {
        let num_controls = control_wires.len() as u8;
        let controled_gate = create_controled_gate(num_controls, &gate.matrix);
        let mut all_wires = vec![0u8;wires.len()+control_wires.len()];
        for i in 0..all_wires.len(){
            if i< wires.len() {
                all_wires[i]=wires[i];
            }else{
                all_wires[i]=control_wires[i-wires.len()];
            }
        }
        let matrix = apply_matrix_to(self.total_qbits, &controled_gate, &all_wires);
        QCS {
            total_qbits: self.total_qbits,
            matrix: matrix.mul(&self.matrix).unwrap(),
        }
    }

    pub fn hadamard(&self, wire: u8) -> QCS {
        let hp = Complex::new(1f64 / 2f64.sqrt(), 0f64);
        self.gate(
            &QCS{total_qbits:1, matrix: Matrix::new(ComplexField, vec![vec![hp, hp], vec![hp, -hp]])},
            &[wire],
        )
    }

    pub fn not(&self, wire: u8) -> QCS {
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        let not_gate = QCS{total_qbits:1, matrix: Matrix::new(ComplexField, vec![vec![z, o],vec![o, z]])};
        self.gate(&not_gate, &[wire])
    }

    pub fn cnot(&self, control: u8, wire: u8) -> QCS {
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        let cnot_gate = QCS{total_qbits:2, matrix: Matrix::new(
            ComplexField,
            vec![
                vec![o, z, z, z],
                vec![z, o, z, z],
                vec![z, z, z, o],
                vec![z, z, o, z],
            ],
        )};
        self.gate(&cnot_gate, &[wire, control])
    }

    pub fn phase(&self, angle:f64) -> QCS {
        let s = Complex::new(angle.cos(), angle.sin());
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        QCS::from(Matrix::new(ComplexField,
            vec![
                vec![o,z],
                vec![z,s]
            ]
        ))
    }

    pub fn swap(&self, x: u8, y: u8) -> QCS {
        let matrix = create_swap_matrix(self.total_qbits, x, y, ComplexField);
        QCS {
            total_qbits: self.total_qbits,
            matrix,
        }
    }

    pub fn run(&self, input: &[Complex<f64>]) -> Vec<Complex<f64>> {
        let mut data = Vec::new();
        for i in 0..input.len() {
            data.push(vec![input[i]]);
        }
        let input_mat = Matrix::new(ComplexField, data);
        let output_mat = self.matrix.mul(&input_mat).unwrap();
        let mut out = Vec::new();
        for i in 0..input.len() {
            out.push(output_mat.value_at(i, 0));
        }
        out
    }

    pub fn value_to_bit_vector(&self, value: usize) -> Vec<u8> {
        let mut v = value;
        let mut result = Vec::new();
        for _ in 0..self.total_qbits {
            result.push((v % 2) as u8);
            v = v / 2;
        }
        result
    }

    pub fn value_to_tensor(&self, value: usize) -> Vec<Complex<f64>> {
        let mut result = Vec::new();
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        for i in 0..pow2(self.total_qbits) {
            if value == i as usize {
                result.push(o);
            } else {
                result.push(z);
            }
        }
        result
    }

    pub fn run_once(&self, input: &[Complex<f64>]) -> Vec<u8> {
        let output = self.run(input);
        let mut rng = rand::thread_rng();
        let probabilities: Vec<f64> = output.into_iter().map(|x| ((x) * (x.conj())).re).collect();

        let y: f64 = rng.gen();
        let mut sum = 0f64;
        for i in 0..input.len() {
            sum += probabilities[i];
            if sum >= y {
                return self.value_to_bit_vector(i);
            }
        }
        return self.value_to_bit_vector(input.len());
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_cnot() {
        let s_not = QCS::new(3).cnot(1, 2).cnot(2, 1).cnot(1, 2);
        let swap = QCS::new(3).swap(2, 1);
        assert_eq!(s_not, swap);
    }

    #[test]
    fn test_controled_gate() {
        let s_not = QCS::new(5).cnot(1, 4).cnot(2, 4);
        let c_not = QCS::new(5).controled_gate(&QCS::new(1).not(0), &[4],&[1])
        .controled_gate(&QCS::new(1).not(0), &[4],&[2]);
        assert_eq!(s_not, c_not);
    }

    #[test]
    fn test_swap() {
        let swap = QCS::new(4).swap(2, 1);
        let res = swap.run_once(&swap.value_to_tensor(0b010));
        let exp = swap.value_to_bit_vector(0b100);
        assert_eq!(res, exp);
    }
    #[test]
    fn test_value_to_tensor() {
        let test = QCS::new(2);
        let res = test.value_to_tensor(0b00);
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        assert_eq!(res, vec![o, z, z, z]);
        let res = test.value_to_tensor(0b01);
        assert_eq!(res, vec![z, o, z, z]);
    }
    #[test]
    fn test_deutsch_0() {
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        let uf = Matrix::new(
            ComplexField,
            vec![
                vec![o, z, z, z],
                vec![z, o, z, z],
                vec![z, z, o, z],
                vec![z, z, z, o],
            ],
        );
        let test = QCS::new(2)
            .hadamard(0)
            .hadamard(1)
            .gate(&QCS::from(uf), &vec![0, 1])
            .hadamard(0)
            .hadamard(1);
        let res = test.run_once(&test.value_to_tensor(0b01));
        assert_eq!(res, vec![1, 0]);
    }

    #[test]
    fn test_deutsch_1() {
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        let uf = Matrix::new(
            ComplexField,
            vec![
                vec![z, o, z, z],
                vec![o, z, z, z],
                vec![z, z, z, o],
                vec![z, z, o, z],
            ],
        );
        let test = QCS::new(2)
            .hadamard(0)
            .hadamard(1)
            .gate(&QCS::from(uf), &vec![0, 1])
            .hadamard(0)
            .hadamard(1);
        let res = test.run_once(&test.value_to_tensor(0b01));
        assert_eq!(res, vec![1, 0]);
    }

    #[test]
    fn test_deutsch_i() {
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        let uf = Matrix::new(
            ComplexField,
            vec![
                vec![o, z, z, z],
                vec![z, o, z, z],
                vec![z, z, z, o],
                vec![z, z, o, z],
            ],
        );
        let test = QCS::new(2)
            .hadamard(0)
            .hadamard(1)
            .gate(&QCS::from(uf), &vec![0, 1])
            .hadamard(0)
            .hadamard(1);
        let res = test.run_once(&test.value_to_tensor(0b01));
        assert_eq!(res, vec![1, 1]);
    }
    #[test]
    fn test_deutsch_not() {
        let z = Complex::new(0f64, 0f64);
        let o = Complex::new(1f64, 0f64);
        let uf = Matrix::new(
            ComplexField,
            vec![
                vec![z, o, z, z],
                vec![o, z, z, z],
                vec![z, z, o, z],
                vec![z, z, z, o],
            ],
        );
        let test = QCS::new(2)
            .hadamard(0)
            .hadamard(1)
            .gate(&QCS::from(uf), &vec![0, 1])
            .hadamard(0)
            .hadamard(1);
        let res = test.run_once(&test.value_to_tensor(0b01));
        assert_eq!(res, vec![1, 1]);
    }
}
