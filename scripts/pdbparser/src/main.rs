#[macro_use]
extern crate clap;

use clap::{App, Arg};
use std::io;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use std::vec::Vec;
use std::collections::HashMap;

fn main() {

    let atomic_radii: Vec<(&str, f32)> = vec![("Ac", 2.15), ("Ag", 1.45), ("Al", 1.21), ("Am", 1.8),
                                            ("Ar", 1.06), ("As", 1.19), ("At", 1.50), ("Au", 1.36),
                                            ("B", 0.84), ("Ba", 2.15), ("Be", 0.96), ("Bh", 1.0),
                                            ("Bi", 1.48), ("Bk", 1.0), ("Br", 1.2), ("C", 0.76),
                                            ("Ca", 1.76), ("Cd", 1.44), ("Ce", 2.04), ("Cf", 1.0),
                                            ("Cl", 1.02), ("Cm", 1.69), ("Co", 1.26), ("Cn", 1.0),
                                            ("Cr", 1.39), ("Cs", 2.44), ("Cu", 1.32), ("Db", 1.0),
                                            ("Ds", 1.0), ("Er", 1.89), ("Es", 1.0), ("Eu", 1.98),
                                            ("F", 0.57), ("Fe", 1.32), ("Fm", 1.0), ("Fr", 2.6),
                                            ("Ga", 1.22), ("Gd", 1.96), ("Ge", 1.20), ("H", 0.31),
                                            ("He", 0.28), ("Hf", 1.75), ("Hg", 1.32), ("Ho", 1.92),
                                            ("Hs", 1.0), ("I", 1.39), ("In", 1.42), ("Ir", 1.41),
                                            ("K", 2.03), ("Kr", 1.16), ("La", 2.07), ("Li", 1.28),
                                            ("Lr", 1.0), ("Lu", 1.87), ("Md", 1.0), ("Mg", 1.41),
                                            ("Mn", 1.39), ("Mo", 1.54), ("Mt", 1.0), ("N", 0.71),
                                            ("Na", 1.66), ("Nb", 1.64), ("Nd", 2.01), ("Ne", 0.58),
                                            ("Ni", 1.24), ("No", 1.0), ("Np", 1.9), ("O", 0.66),
                                            ("Os", 1.44), ("P", 1.07), ("Pa", 2.0), ("Pb", 1.46),
                                            ("Pd", 1.39), ("Pm", 1.99), ("Po", 1.40), ("Pr", 2.03),
                                            ("Pt", 1.36), ("Pu", 1.87), ("Ra", 2.21), ("Rb", 2.2),
                                            ("Re", 1.51), ("Rf", 1.0), ("Rg", 1.0), ("Rh", 1.42),
                                            ("Rn", 1.0), ("Ru", 1.46), ("S", 1.05), ("Sb", 1.39),
                                            ("Sc", 1.7), ("Se", 1.2), ("Sg", 1.0), ("Si", 1.11),
                                            ("Sm", 1.98), ("Sn", 1.39), ("Sr", 1.95), ("Ta", 1.7),
                                            ("Tb", 1.94), ("Tc", 1.47), ("Te", 1.38), ("Th", 2.06),
                                            ("Ti", 1.6), ("Tl", 1.45), ("Tm", 1.90), ("U", 1.96),
                                            ("V", 1.53), ("W", 1.62), ("Xe", 1.40), ("Y", 1.9),
                                            ("Yb", 1.87), ("Zn", 1.22), ("Zr", 1.75)];

    let atomic_radii: HashMap<_, _> = atomic_radii.into_iter().collect();

    //Read in the filename and the scaling required.
    let matches = App::new("PDB Parser")
                      .about("Converts a given pdb file into an edge list, with atoms as nodes and
                      a specified scaling")
                      .arg(Arg::with_name("filename")
                           .help("the location of the PDB file")
                           .required(true))
                      .arg(Arg::with_name("scaling")
                           .short("s")
                           .long("scaling")
                           .takes_value(true)
                           .help("Sets the cutoff radius for edge generation"))
                      .get_matches();

      let filename: String = matches.value_of("filename").unwrap().to_string();
      let scaling: f32 = value_t!(matches, "scaling", f32).unwrap_or(2.5);

      // Split the string, removing the ".pdb" extension and creating a .dat file
      let output = format!("{}.{}.dat", filename[0..filename.len()-4].to_string(), scaling);
      //println!("Generating edgelist: {}",output);

      // Read the file contents
      let path = Path::new(&filename);
      let display = path.display();
      let f = match File::open(&path) {
          Ok(file) => file,
          Err(e) => panic!("couldn't read {}: {}", display, e)
      };
      let file = io::BufReader::new(&f);
      let mut sizes: Vec<(Vec<f32>, String)> = Vec::new();
      for line in file.lines() {
          let l = line.unwrap();
          // if record is an atom, push the x,y,z co-ordinates, and the element, to a file.
          if l[0..4].to_string() == "ATOM" {
              let record: Vec<f32>= vec![l[30..38].trim().parse::<f32>().unwrap(), l[38..46].trim().parse::<f32>().unwrap(), l[46..54].trim().parse::<f32>().unwrap()];

              let element = l[76..78].trim().to_string();
              sizes.push((record,element));
          }
      }

      // Now we have the list of atomic co-ordinates; for every atom, get the distance between its
      // neighbours and output the weighted edge.

      let output = Path::new(&output);
      let display = output.display();
      let f = match File::create(&output) {
          Ok(file) => file,
          Err(e) => panic!("couldn't create {}: {}", display, e)
      };
      let mut writer = io::BufWriter::new(&f);

      for (i, row) in sizes.iter().enumerate() {
          for (j,row2) in sizes[0..i].iter().enumerate() {
              let distance_squared: f32 = get_distance_squared(&row.0, &row2.0);
              let elem: &str= &row.1[..];
              let elem2: &str = &row2.1[..];
              let cutoff = (atomic_radii.get(elem).unwrap() + atomic_radii.get(elem2).unwrap()) * scaling;
              if distance_squared < cutoff*cutoff {
                  let distance = distance_squared.sqrt();
                  let weight = (cutoff - distance) / cutoff;
                  write!(&mut writer, "{} {} {}\n", i+1, j+1, weight);
              }
          }
      }
}

fn get_distance_squared(row:&Vec<f32>, row2: &Vec<f32>) -> f32 {

    (row[0]-row2[0]).powi(2) + (row[1]-row2[1]).powi(2) + (row[2]-row2[2]).powi(2)

}
