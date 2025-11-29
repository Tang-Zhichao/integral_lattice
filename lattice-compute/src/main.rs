use std::io::{self, Write};
use std::time::Instant;

struct Matrix {
    data: Vec<Vec<i64>>,
    n: usize,
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let temp = b;
        b = a % b;
        a = temp;
    }
    a
}

fn extended_gcd(a: i64, b: i64) -> (i64, i64, i64) {
    if b == 0 {
        return (a, 1, 0);
    }
    let (g, x1, y1) = extended_gcd(b, a % b);
    (g, y1, x1 - (a / b) * y1)
}

impl Matrix {
    fn new(n: usize) -> Self {
        Matrix {
            data: vec![vec![0; n]; n],
            n,
        }
    }

    fn input(&mut self) -> io::Result<()> {
        println!(
            "Enter {}×{} matrix row by row (space-separated integers):",
            self.n, self.n
        );

        for i in 0..self.n {
            print!("Row {}: ", i + 1);
            io::stdout().flush()?;

            let mut line = String::new();
            io::stdin().read_line(&mut line)?;

            let values: Result<Vec<i64>, _> =
                line.trim().split_whitespace().map(|s| s.parse()).collect();

            match values {
                Ok(row) if row.len() == self.n => {
                    self.data[i] = row;
                }
                Ok(_) => {
                    println!("Error: Expected {} values", self.n);
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "Wrong number of values",
                    ));
                }
                Err(e) => {
                    println!("Error parsing input: {}", e);
                    return Err(io::Error::new(io::ErrorKind::InvalidInput, "Parse error"));
                }
            }
        }

        Ok(())
    }

    fn print(&self) {
        for row in &self.data {
            for val in row {
                print!("{:6} ", val);
            }
            println!();
        }
    }

    fn is_symmetric(&self) -> bool {
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                if self.data[i][j] != self.data[j][i] {
                    return false;
                }
            }
        }
        true
    }

    fn determinant(&self) -> i64 {
        if self.n == 1 {
            return self.data[0][0];
        }
        if self.n == 2 {
            return self.data[0][0] * self.data[1][1] - self.data[0][1] * self.data[1][0];
        }

        let mut temp: Vec<Vec<f64>> = self
            .data
            .iter()
            .map(|row| row.iter().map(|&x| x as f64).collect())
            .collect();

        let mut det = 1.0;

        for i in 0..self.n {
            let mut pivot = i;
            for j in (i + 1)..self.n {
                if temp[j][i].abs() > temp[pivot][i].abs() {
                    pivot = j;
                }
            }

            if temp[pivot][i].abs() < 1e-10 {
                return 0;
            }

            if pivot != i {
                temp.swap(i, pivot);
                det = -det;
            }

            det *= temp[i][i];

            for j in (i + 1)..self.n {
                let factor = temp[j][i] / temp[i][i];
                for k in i..self.n {
                    temp[j][k] -= temp[i][k] * factor;
                }
            }
        }

        det.round() as i64
    }

    // Compute Smith Normal Form to find invariant factors
    fn smith_normal_form(&self, modulo: i64) -> Vec<i64> {
        let mut mat = self.data.clone();
        let m = self.n;
        let n = self.n;

        // Apply modulo to all entries
        for i in 0..m {
            for j in 0..n {
                mat[i][j] = ((mat[i][j] % modulo) + modulo) % modulo;
            }
        }

        let mut rank = 0;

        for col in 0..n.min(m) {
            // Find non-zero pivot
            let mut pivot_row = None;
            for row in rank..m {
                if mat[row][col] != 0 {
                    pivot_row = Some(row);
                    break;
                }
            }

            if pivot_row.is_none() {
                continue;
            }

            let pivot_row = pivot_row.unwrap();
            if pivot_row != rank {
                mat.swap(rank, pivot_row);
            }

            // Reduce all other entries in this column
            for row in 0..m {
                if row == rank {
                    continue;
                }

                if mat[row][col] == 0 {
                    continue;
                }

                let g = gcd(mat[rank][col], mat[row][col]);
                if g == 0 {
                    continue;
                }

                let (_, s, t) = extended_gcd(mat[rank][col], mat[row][col]);
                let a = mat[rank][col] / g;
                let b = mat[row][col] / g;

                for col_idx in 0..n {
                    let temp1 = mat[rank][col_idx];
                    let temp2 = mat[row][col_idx];
                    mat[rank][col_idx] = ((s * temp1 + t * temp2) % modulo + modulo) % modulo;
                    mat[row][col_idx] = ((-b * temp1 + a * temp2) % modulo + modulo) % modulo;
                }
            }

            rank += 1;
        }

        // Extract diagonal elements
        let mut invariants = Vec::new();
        for i in 0..rank.min(n) {
            let val = mat[i][i] % modulo;
            if val > 0 {
                invariants.push(gcd(val, modulo));
            }
        }

        // Ensure divisibility chain
        invariants.sort();
        let mut cleaned = Vec::new();
        for &inv in &invariants {
            if inv > 1 {
                cleaned.push(inv);
            }
        }

        cleaned
    }

    fn multiply_mod(&self, vec: &[i64], modulo: i64) -> Vec<i64> {
        let mut result = vec![0i64; self.n];
        for i in 0..self.n {
            let mut sum: i128 = 0;
            for j in 0..self.n {
                sum += (self.data[i][j] as i128) * (vec[j] as i128);
            }
            result[i] = ((sum % modulo as i128 + modulo as i128) % modulo as i128) as i64;
        }
        result
    }

    fn find_modular_null_vectors(&self, modulo: i64, show_all: bool, display_time_limit: f64) {
        println!(
            "\nFinding all vectors v where (A @ v) ≡ 0 (mod {})...",
            modulo
        );
        if show_all {
            println!(
                "(Will display vectors for {:.1} seconds, then continue counting)",
                display_time_limit
            );
        }

        let start = Instant::now();
        let mut vec = vec![0i64; self.n];
        let mut count = 0;
        let mut solutions = Vec::new();
        let mut stopped_showing = false;

        self.find_null_vectors_recursive(
            &mut vec,
            0,
            modulo,
            &mut count,
            &mut solutions,
            show_all,
            &start,
            display_time_limit,
            &mut stopped_showing,
        );

        let duration = start.elapsed();

        println!("\n{}", "=".repeat(70));
        println!("RESULTS:");
        println!("Total vectors found: {}", count);
        println!("Expected (|det|): {}", modulo.abs());
        println!(
            "Match: {}",
            if count == modulo.abs() as usize {
                "✓ YES"
            } else {
                "✗ NO"
            }
        );
        println!("Time elapsed: {:.6} seconds", duration.as_secs_f64());

        if count > 0 {
            let rate = count as f64 / duration.as_secs_f64();
            println!("Search rate: {:.0} vectors/second", rate);
        }

        println!("{}", "=".repeat(70));
    }

    fn find_null_vectors_recursive(
        &self,
        vec: &mut Vec<i64>,
        pos: usize,
        modulo: i64,
        count: &mut usize,
        solutions: &mut Vec<Vec<i64>>,
        show_all: bool,
        start_time: &Instant,
        display_time_limit: f64,
        stopped_showing: &mut bool,
    ) {
        if pos == self.n {
            let result = self.multiply_mod(vec, modulo);

            if result.iter().all(|&x| x == 0) {
                *count += 1;

                if show_all {
                    let elapsed = start_time.elapsed().as_secs_f64();

                    if elapsed < display_time_limit {
                        print!("{}. [", count);
                        for (i, &val) in vec.iter().enumerate() {
                            print!("{}", val);
                            if i < self.n - 1 {
                                print!(", ");
                            }
                        }
                        println!("]");
                    } else if !*stopped_showing {
                        println!(
                            "\n... ({:.1}s reached, continuing search without display) ...\n",
                            display_time_limit
                        );
                        *stopped_showing = true;
                    }
                }

                if solutions.len() < 1000 {
                    solutions.push(vec.clone());
                }
            }
            return;
        }

        for val in 0..modulo {
            vec[pos] = val;
            self.find_null_vectors_recursive(
                vec,
                pos + 1,
                modulo,
                count,
                solutions,
                show_all,
                start_time,
                display_time_limit,
                stopped_showing,
            );
        }
    }
}

fn main() -> io::Result<()> {
    println!("{}", "=".repeat(70));
    println!("Matrix Modular Null Space Finder with Group Structure Analysis");
    println!("{}", "=".repeat(70));

    print!("\nEnter matrix dimension (n for n×n matrix): ");
    io::stdout().flush()?;

    let mut input = String::new();
    io::stdin().read_line(&mut input)?;

    let n: usize = match input.trim().parse() {
        Ok(num) if num > 0 && num <= 20 => num,
        _ => {
            println!("Invalid dimension. Please use 1 ≤ n ≤ 20");
            return Ok(());
        }
    };

    let mut matrix = Matrix::new(n);
    matrix.input()?;

    println!("\nInput Matrix:");
    matrix.print();

    let check_start = Instant::now();
    let symmetric = matrix.is_symmetric();
    println!("\nIs symmetric: {}", if symmetric { "Yes" } else { "No" });

    let det = matrix.determinant();
    println!("Determinant: {}", det);
    let check_time = check_start.elapsed();
    println!("(Check time: {:.6} seconds)", check_time.as_secs_f64());

    if symmetric && det != 0 {
        let modulo = det.abs();
        println!("\nModulo: {}", modulo);

        // Compute Smith Normal Form / Invariant Factors
        println!("\n{}", "=".repeat(70));
        println!("GROUP STRUCTURE ANALYSIS:");
        println!("{}", "=".repeat(70));

        let snf_start = Instant::now();
        let invariants = matrix.smith_normal_form(modulo);
        let snf_time = snf_start.elapsed();

        println!("Computing invariant factors (Smith Normal Form)...");
        println!("Time: {:.6} seconds", snf_time.as_secs_f64());

        if invariants.is_empty() {
            println!(
                "\nInvariant factors: [{}] (free group Z/{}Z)",
                modulo, modulo
            );
            println!("Group structure: Z/{} ≅ Z/{}", modulo, modulo);
        } else {
            println!("\nInvariant factors: {:?}", invariants);
            print!("Group structure: ");
            for (i, &inv) in invariants.iter().enumerate() {
                print!("Z/{}", inv);
                if i < invariants.len() - 1 {
                    print!(" ⊕ ");
                }
            }
            println!();

            let product: i64 = invariants.iter().product();
            println!("Product of invariants: {}", product);
            println!("Verification: product should equal |det| = {}", modulo);
        }

        let estimated_ops = modulo.pow(n as u32);
        if estimated_ops > 100_000_000 {
            print!(
                "\nWarning: Will check ~{} combinations. Continue? (y/n): ",
                estimated_ops
            );
            io::stdout().flush()?;

            let mut choice = String::new();
            io::stdin().read_line(&mut choice)?;

            if !choice.trim().eq_ignore_ascii_case("y") {
                return Ok(());
            }
        }

        print!("\nShow solutions? (y/n): ");
        io::stdout().flush()?;
        let mut choice = String::new();
        io::stdin().read_line(&mut choice)?;
        let show_all = choice.trim().eq_ignore_ascii_case("y");

        let display_time_limit = if show_all {
            print!("Display vectors for how many seconds? (e.g., 5, 30, 60): ");
            io::stdout().flush()?;
            let mut time_input = String::new();
            io::stdin().read_line(&mut time_input)?;
            time_input.trim().parse().unwrap_or(10.0)
        } else {
            0.0
        };

        matrix.find_modular_null_vectors(modulo, show_all, display_time_limit);
    } else if !symmetric {
        println!("\nMatrix is not symmetric. No further computation performed.");
    } else {
        println!("\nMatrix is degenerate (determinant = 0). No further computation performed.");
    }

    Ok(())
}
