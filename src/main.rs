extern crate linearalgebra;
use crate::linearalgebra::Ring;
use linearalgebra::Matrix;
use linearalgebra::field::ComplexField;
use num_complex::Complex;
use std::num;

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

fn linear_search(list:&[u8], value:u8) -> Option<u8> {
    for i in 0..list.len() {
        if list[i] == value {
            return Some(i as u8)
        }
    }
    return None;
}

//returns a list of swaps required to send the entries to the position 0,1,..
fn create_swap_sequence(total_qbits:u8, entries:&[u8], ) -> Vec<(u8,u8)> {
    let mut swaps = Vec::new();
    let mut currect_pos = Vec::new();
    for i in 0..total_qbits {
        currect_pos.push(i);
    }
    let mut next = 0;
    for e in entries {
        if let Some(pos) = linear_search(&currect_pos[..], *e) {
            if next != pos {
                swaps.push((pos, next));
                currect_pos.swap(pos as usize, next as usize);
            }

            next+=1;
        }
    }

    swaps
}

pub fn  create_swap_matrix(total_qbits:u8, qbit_x:u8,qbit_y:u8, complex_field: ComplexField) -> Matrix<ComplexField>{
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
    let matrix = Matrix::new(complex_field, mat);
    matrix.transpose()
}

pub fn create_swap_matrix_for_entries(total_qbits:u8, entries:&[u8]) -> Matrix<ComplexField> {
    let complex_field = ComplexField;
    let swap_seq = create_swap_sequence(total_qbits, entries);
    let mut prod = create_swap_matrix(total_qbits,0,0,complex_field);

    for swap in swap_seq {
        prod = create_swap_matrix(total_qbits,swap.0,swap.1,complex_field).mul(&prod).unwrap();
    }
    prod
}
fn stack_mat(total_qbits:u8, matrix: Matrix<ComplexField>) -> Matrix<ComplexField> {
    let dim = pow2(total_qbits);
    if dim == matrix.rows() as u8{
        matrix
    }else {
        stack_mat(total_qbits, stack_twice(matrix))
    }
}
fn stack_twice(matrix: Matrix<ComplexField>) -> Matrix<ComplexField> {
    let ring = ComplexField;
    let rows = 2*matrix.rows();
    let mut data = vec![vec![ComplexField.zero();rows];rows];
    for i in 0..rows{
        for j in 0..rows{
            if i< rows/2 {
                if j<rows/2 {
                    data[i][j] = matrix.value_at(i, j);
                }
            } else{
                if j<rows/2 {

                } else{
                    data[i][j] = matrix.value_at(i-rows/2, j-rows/2);
                }
            }
        }
    }
    Matrix::new(ring, data)
}
pub fn conjugate_transpose(matrix: &Matrix<ComplexField>) -> Matrix<ComplexField> {
    let ring = ComplexField;
    let rows = matrix.rows();
    let mut data = vec![vec![ComplexField.zero();rows];rows];
    for i in 0..rows{
        for j in 0..rows{
            data[i][j]=matrix.value_at(j, i).conj();
        }
    }
    Matrix::new(ring, data)
}
pub fn apply_matrix_to(total_qbits:u8, matrix: &Matrix<ComplexField>, bits_to_apply2:&[u8]) -> Matrix<ComplexField> {
    let swaps = create_swap_matrix_for_entries(total_qbits, bits_to_apply2);
    let inv_swaps = conjugate_transpose(&swaps);
    let ramped_mat = stack_mat(total_qbits, matrix.clone());
    inv_swaps.mul(&ramped_mat.mul(&swaps).unwrap()).unwrap()
}
fn print_real_matrix(mat:&Matrix<ComplexField>) {
    for i in 0..mat.rows() {
        for j in 0..mat.columns() {
            print!("{:.2?} ", mat.value_at(i, j).re);
        }
        println!();
    }
}

fn main() {
    let ring = ComplexField;
    let total_qbits = 4;
    let qbit_x = 1;
    let qbit_y = 2;
    let complex_field = ComplexField;

    let mat = create_swap_matrix(total_qbits, qbit_x, qbit_y, complex_field);

    let swap_seq = create_swap_sequence(total_qbits, &vec![0,1]);
    //print_real_matrix(&mat);
    println!("{:?}", swap_seq);
    let hp = Complex::new(1.0/(2.0 as f64).sqrt(), 0f64);
    let mat = Matrix::new(ring.clone(),vec![vec![hp,hp], vec![hp, -hp]]);
    let had2= apply_matrix_to(2, &mat, &vec![0]);
    let init = Matrix::new(ring.clone(),vec![vec![Complex::new(0f64, 0f64)],
    vec![Complex::new(0f64, 0f64)], vec![Complex::new(0f64, 0f64)], vec![Complex::new(1f64, 0f64)]]);
    print_real_matrix(&had2);
    print_real_matrix(&had2.mul(&init).unwrap());
}
