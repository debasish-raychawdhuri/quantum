extern crate linearalgebra;


fn pow2(x:u8) -> u8 {
    1<<x
}

fn get_column_for_one(size:u8, one_position: u8) -> Vec<u8> {
    let mut v = Vec::new();
    for i in 0..size {
        if i == one_position {
            v.push(1u8);
        }else{
            v.push(0u8);
        }
    }
    v

}

pub fn create_swap_matrix(total_qbits:u8, qbit_x:u8,qbit_y:u8) -> Vec<Vec<u8>>{
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
    mat
}


fn main() {
    let total_qbits = 4;
    let qbit_x = 1;
    let qbit_y = 2;

    let mat = create_swap_matrix(total_qbits, qbit_x, qbit_y);
    for i in 0..mat.len() {
        for j in  0..mat.len() {
            print!("{}",mat[j][i]);
            if j < mat.len()-1{
                print!(",");
            }
        }
        if i < mat.len() -1 {
            print!(";");
        }

    }
    println!("");

}
