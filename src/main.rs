use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Serialize, Deserialize)]
struct Point {
    x: f64,
    z: f64,
}

#[derive(Serialize, Deserialize)]
struct SoilLayer {
    soilcode: String,
    points: Vec<Point>,
}

#[derive(Serialize, Deserialize)]
struct BBFSettings {
    top: f64,
    bottom: f64,
    gridsize_x: f64,
    gridsize_z: f64,
    tangent_top: f64,
    tangent_bottom: f64,
    tangent_gridsize_z: f64,
}

impl SoilLayer {
    fn xmin(&self) -> f64 {
        let mut smallest = 1.0e-9;

        // Iterate over the rest of the points
        for point in &self.points {
            if point.x < smallest {
                smallest = point.x;
            }
        }

        smallest
    }
    fn zmin(&self) -> f64 {
        let mut smallest = 1.0e-9;

        // Iterate over the rest of the points
        for point in &self.points {
            if point.z < smallest {
                smallest = point.z;
            }
        }

        smallest
    }
    fn xmax(&self) -> f64 {
        let mut largest = -1.0e-9;

        // Iterate over the rest of the points
        for point in &self.points {
            if point.x > largest {
                largest = point.x;
            }
        }

        largest
    }
    fn zmax(&self) -> f64 {
        let mut largest = -1.0e-9;

        // Iterate over the rest of the points
        for point in &self.points {
            if point.z > largest {
                largest = point.z;
            }
        }

        largest
    }
}

#[derive(Serialize, Deserialize)]
struct Soil {
    id: i64,
    code: String,
    yd: f64,
    ys: f64,
    c: f64,
    phi: f64,
}

#[derive(Serialize, Deserialize)]
struct Geometry {
    gridsize_x: f64,
    gridsize_z: f64,
    layers: Vec<SoilLayer>,
    phreatic_line: Vec<Point>,
    soils: Vec<Soil>,
    bbf_settings: BBFSettings,
}

impl Geometry {
    fn xmin(&self) -> f64 {
        let mut smallest = 1.0e9;

        for layer in self.layers.iter() {
            smallest = layer.xmin().min(smallest);
        }

        smallest
    }
    fn xmax(&self) -> f64 {
        let mut largest = -1.0e9;

        for layer in self.layers.iter() {
            largest = layer.xmax().max(largest);
        }

        largest
    }
    fn zmin(&self) -> f64 {
        let mut smallest = 1.0e9;

        for layer in self.layers.iter() {
            smallest = layer.zmin().min(smallest);
        }

        smallest
    }
    fn zmax(&self) -> f64 {
        let mut largest = -1.0e9;

        for layer in self.layers.iter() {
            largest = layer.zmax().max(largest);
        }

        largest
    }

    fn soil_id(&self, soilcode: &String) -> i64 {
        for soil in &self.soils {
            if soil.code == *soilcode {
                return soil.id;
            }
        }
        -1
    }

    fn to_soil_matrix(&self) -> Vec<Vec<i64>> {
        // get the width and height of the model
        let dx = self.xmax() - self.xmin();
        let dz = self.zmax() - self.zmin();

        // calculate the number of elements
        let numx = (dx / self.gridsize_x).ceil() as usize;
        let numz = (dz / self.gridsize_z).ceil() as usize;

        // create a matrix for the result
        // the matrix will be x, z aranged
        let mut matrix = vec![vec![0; numz]; numx];

        // fill in the values
        let xmin = self.xmin();
        let zmax = self.zmax();
        // println!(
        //     "xmin: {}, xmax: {}, zmin: {}, zmax: {}, numx: {}, numz: {}",
        //     xmin, xmax, zmin, zmax, numx, numz
        // );
        for i in 0..numx {
            for j in 0..numz {
                let x = xmin + (i as f64 + 0.5) * self.gridsize_x;
                let z = zmax - (j as f64 + 0.5) * self.gridsize_z;
                //println!("x: {}, z: {}", x, z);

                for layer in &self.layers {
                    let p = Point { x: x, z: z };
                    if is_point_in_polygon(&p, &layer.points) {
                        let soil_id = self.soil_id(&layer.soilcode);
                        // println!(
                        //     "i: {}, j: {}, x: {}, z: {}, soil_id: {}",
                        //     i, j, x, z, soil_id
                        // );
                        matrix[i][j] = soil_id;
                        break;
                    }
                }
            }
        }
        matrix
    }

    fn soil_by_id(&self, id: i64) -> Option<&Soil> {
        for soil in &self.soils {
            if soil.id == id {
                return Some(soil);
            }
        }
        None
    }
}

struct Stresses {
    m_stot: Vec<Vec<f64>>,
    m_u: Vec<Vec<f64>>,
    m_seff: Vec<Vec<f64>>,
}

impl Stresses {
    fn new(numx: usize, numz: usize) -> Stresses {
        let m_stot: Vec<Vec<f64>> = vec![vec![0.0; numz]; numx];
        let m_u: Vec<Vec<f64>> = vec![vec![0.0; numz]; numx];
        let m_seff: Vec<Vec<f64>> = vec![vec![0.0; numz]; numx];

        Stresses {
            m_stot,
            m_u,
            m_seff,
        }
    }
}

struct Analysis {
    geometry: Geometry,
    m_soil: Vec<Vec<i64>>,
    m_phreatic: Vec<Vec<bool>>,
    stresses: Stresses,
}

impl Analysis {
    fn new(geometry: Geometry) -> Analysis {
        let m_soil = &geometry.to_soil_matrix();
        let numx = m_soil.len();
        let numz = m_soil[0].len();

        let m_phreatic: Vec<Vec<bool>> = vec![vec![false; numz]; numx];

        Analysis {
            geometry,
            m_soil: m_soil.clone(),
            m_phreatic,
            stresses: Stresses::new(numx, numz),
        }
    }

    fn generate_stresses_matrix(&self) {}

    fn execute_bishop(&mut self) {
        let numx = self.m_soil.len();
        let numz = self.m_soil[0].len();

        let mut m_tau: Vec<Vec<f64>> = vec![vec![0.0; numz]; numx];

        for i in 0..numx {
            let mut s = 0.0;
            let mut u = 0.0;
            for j in 0..numz {
                if self.m_phreatic[i][j] {
                    u += self.geometry.gridsize_z * 10.0;
                }
                // println!(
                //     "i = {}, j = {}, p = {}, u = {}",
                //     i, j, self.m_phreatic[i][j], u
                // );

                if self.m_soil[i][j] != 0 {
                    let soil = self.geometry.soil_by_id(self.m_soil[i][j]).unwrap();
                    if self.m_phreatic[i][j] {
                        s += self.geometry.gridsize_z * soil.ys;
                    } else {
                        s += self.geometry.gridsize_z * soil.yd;
                    }

                    let seff = (s - u).max(0.0);
                    m_tau[i][j] = soil.c + seff * soil.phi.to_radians().tan();
                }

                self.stresses.m_u[i][j] = u;
                self.stresses.m_stot[i][j] = s;
                self.stresses.m_seff[i][j] = (s - u).max(0.0);
            }
        }
        print_matrix(&m_tau);
    }
}

fn is_point_in_polygon(p: &Point, polygon: &[Point]) -> bool {
    let mut crossings = 0;
    let n = polygon.len();

    for i in 0..n {
        let j = (i + 1) % n;
        let (x1, z1) = (polygon[i].x, polygon[i].z);
        let (x2, z2) = (polygon[j].x, polygon[j].z);

        if (p.z >= z1 && p.z < z2 || p.z >= z2 && p.z < z1)
            && (p.x < (x2 - x1) * (p.z - z1) / (z2 - z1) + x1)
        {
            crossings += 1;
        } else if p.x == x1 && p.z == z1 {
            // Point is exactly on a vertex
            return true;
        }
    }

    crossings % 2 != 0
}

fn print_matrix(matrix: &Vec<Vec<f64>>) {
    let rows = matrix.len();
    let cols = matrix[0].len();

    for j in 0..cols {
        for i in (0..rows).rev() {
            print!("{:6.1} ", matrix[i][j]); // Adjust formatting as needed
        }
        println!(); // Move to the next line for the next column
    }
}

fn main() {
    let contents = fs::read_to_string("testdata/geometry.json").expect("Unable to read file");

    let geometry: Geometry = serde_json::from_str(&contents).expect("Failed to parse JSON");
    let mut analysis = Analysis::new(geometry);

    analysis.execute_bishop();

    print_matrix(&analysis.stresses.m_stot);

    //println!("M[0][0]: {}", M[0][0])
}
