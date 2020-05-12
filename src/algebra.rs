extern crate linearalgebra;
use linearalgebra::field::ComplexField;
use linearalgebra::Matrix;
use linearalgebra::Ring;
use num_complex::Complex;
use std::num;

pub fn pow2(x: u8) -> u64 {
    1u64 << (x as u64)
}

fn get_column_for_one(size: u64, one_position: u64) -> Vec<Complex<f64>> {
    let mut v = Vec::new();
    for i in 0..size {
        if i == one_position {
            v.push(Complex::new(1f64, 0f64));
        } else {
            v.push(Complex::new(0f64, 0f64));
        }
    }
    v
}

fn linear_search(list: &[u8], value: u8) -> Option<u8> {
    for i in 0..list.len() {
        if list[i] == value {
            return Some(i as u8);
        }
    }
    return None;
}

//returns a list of swaps required to send the entries to the position 0,1,..
fn create_swap_sequence(total_qbits: u8, entries: &[u8]) -> Vec<(u8, u8)> {
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

            next += 1;
        }
    }

    swaps
}

pub fn create_swap_matrix(
    total_qbits: u8,
    qbit_x: u8,
    qbit_y: u8,
    complex_field: ComplexField,
) -> Matrix<ComplexField> {
    let dim = pow2(total_qbits);
    let x_rep = 1 << qbit_x;
    let y_rep = 1 << qbit_y;
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

pub fn create_swap_matrix_for_entries(total_qbits: u8, entries: &[u8]) -> Matrix<ComplexField> {
    let complex_field = ComplexField;
    let swap_seq = create_swap_sequence(total_qbits, entries);
    let mut prod = create_swap_matrix(total_qbits, 0, 0, complex_field);

    for swap in swap_seq {
        prod = create_swap_matrix(total_qbits, swap.0, swap.1, complex_field)
            .mul(&prod)
            .unwrap();
    }
    prod
}

fn stack_mat(total_qbits: u8, matrix: Matrix<ComplexField>) -> Matrix<ComplexField> {
    if matrix.rows() != matrix.columns() {
        panic!("Only square matrix supported");
    }
    let ring = ComplexField;
    let rows = pow2(total_qbits) as usize;
    let mut data = vec![vec![ComplexField.zero(); rows]; rows];
    for i in 0..rows {
        let row_block = i / matrix.rows();
        for j in 0..rows {
            let col_block = j / matrix.rows();

            if row_block == col_block {
                data[i][j] =
                    matrix.value_at(i - matrix.rows() * col_block, j - matrix.rows() * col_block);
            }
        }
    }
    Matrix::new(ring, data)
}

fn stack_id(total_qbits: u8, matrix: Matrix<ComplexField>) -> Matrix<ComplexField> {
    if matrix.rows() != matrix.columns() {
        panic!("Only square matrix supported");
    }
    let ring = ComplexField;
    let rows = pow2(total_qbits) as usize;
    let total_blocks = rows / matrix.rows();
    let mut data = vec![vec![ComplexField.zero(); rows]; rows];
    for i in 0..rows {
        let row_block = i / matrix.rows();
        for j in 0..rows {
            let col_block = j / matrix.rows();
            if row_block == col_block && col_block == total_blocks - 1 {
                data[i][j] =
                    matrix.value_at(i - matrix.rows() * col_block, j - matrix.rows() * col_block);
            } else if i == j {
                data[i][j] = ring.one();
            }
        }
    }
    Matrix::new(ring, data)
}

pub fn conjugate_transpose(matrix: &Matrix<ComplexField>) -> Matrix<ComplexField> {
    let ring = ComplexField;
    let rows = matrix.rows();
    let mut data = vec![vec![ComplexField.zero(); rows]; rows];
    for i in 0..rows {
        for j in 0..rows {
            data[i][j] = matrix.value_at(j, i).conj();
        }
    }
    Matrix::new(ring, data)
}
pub fn apply_matrix_to(
    total_qbits: u8,
    matrix: &Matrix<ComplexField>,
    bits_to_apply2: &[u8],
) -> Matrix<ComplexField> {
    let swaps = create_swap_matrix_for_entries(total_qbits, bits_to_apply2);
    let inv_swaps = conjugate_transpose(&swaps);
    let ramped_mat = stack_mat(total_qbits, matrix.clone());
    inv_swaps.mul(&ramped_mat.mul(&swaps).unwrap()).unwrap()
}
fn print_real_matrix(mat: &Matrix<ComplexField>) {
    for i in 0..mat.rows() {
        for j in 0..mat.columns() {
            print!("{:.2?} ", mat.value_at(i, j).re);
        }
        println!();
    }
}
