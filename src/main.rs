extern crate linearalgebra;
use linearalgebra::Matrix;
use linearalgebra::field::ComplexField;
use num_complex::Complex;

fn pow2(x:u8) -> u8 {
    1<<x
}

fn get_column_for_one(size:u8, one_position: u8) -> Vec<Complex<f64>> {
    let mut v = Vec::new();
    for i in 0..size {
        if i == one_position {
            v.push(Complex::new(1f64, 0f64));
        }else{
            v.push(Complex::new(0f64, 0f64));
        }
    }
    v

}

pub fn  create_swap_matrix(total_qbits:u8, qbit_x:u8,qbit_y:u8, complexField: ComplexField) -> Matrix<ComplexField>{
    let dim = pow2(total_qbits);
    let x_rep = 1<<qbit_x;
    let y_rep = 1<<qbit_y;
    let mut mat = Vec::new();
    for i in 0..dim {
        let bit_x = (i & x_rep) >> qbit_x;
        let bit_y = (i & y_rep) >> qbit_y;
        let to_x = bit_x << qbit_y;
        let to_y = bit_y << qbit_x;
        let shift_x = (i ^ (i & x_rep)) | to_y;
        let shift_y = (shift_x ^ (i & y_rep)) | to_x;
        mat.push(get_column_for_one(dim, shift_y));
    }
    let matrix = Matrix::new(complexField, mat);
    matrix.transpose()
}

fn print_real_matrix(mat:&Matrix<ComplexField>) {
    for i in 0..mat.rows() {
        for j in 0..mat.columns() {
            print!("{:?} ", mat.value_at(i, j).re);
        }
        println!();
    }
}

fn main() {
    let total_qbits = 4;
    let qbit_x = 1;
    let qbit_y = 2;
    let complexField = ComplexField;

    let mat = create_swap_matrix(total_qbits, qbit_x, qbit_y, complexField);

    print_real_matrix(&mat);

}
