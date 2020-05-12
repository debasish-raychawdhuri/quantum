
use quantum_bits::*;
use linearalgebra::field::ComplexField;
use linearalgebra::Matrix;
use num_complex::Complex;

fn main() {
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
    println!("The result is {:?}", res);
}
