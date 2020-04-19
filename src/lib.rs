mod algebra;
use linearalgebra::Matrix;
use linearalgebra::field::ComplexField;
use algebra::*;
use num_complex::Complex;

#[derive(PartialEq, Clone, Debug)]
pub struct QCS{
    total_qbits:u8,
    matrix: Matrix<ComplexField>
}


impl QCS {
    pub fn new(total_qbits:u8) -> Self {
        let complex_field = ComplexField;
        let matrix = create_swap_matrix(total_qbits, 0, 0, complex_field);
        QCS {
            total_qbits,
            matrix
        }
    }

    pub fn gate(&self, gate:Matrix<ComplexField>, wires: &[u8]) -> QCS {
        let matrix = apply_matrix_to(self.total_qbits, &gate, wires);
        QCS {
            total_qbits:self.total_qbits,
            matrix: matrix.mul(&self.matrix).unwrap()
        }
    }

    pub fn hadamard(&self, wire:u8) -> QCS {
        let hp = Complex::new(1f64/2f64.sqrt(), 0f64);
        self.gate(Matrix::new(ComplexField, vec![vec![hp,hp],vec![hp,-hp]]),&[wire])
    }

    pub fn not(&self, wire:u8) -> QCS {
        let z = Complex::new(0f64,0f64);
        let o = Complex::new(1f64, 0f64);
        let not_gate = Matrix::new(ComplexField,vec![vec![z,o,o,z]]);
        self.gate(not_gate,&[wire])
    }
}
