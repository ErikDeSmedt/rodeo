#![deny(non_camel_case_types)]

extern crate nalgebra as na;
extern crate alga;

mod methods;
//mod problems;

use na::{Matrix6,Vector6};

use alga::linear::NormedSpace;
use alga::general::RealField;

use crate::methods::ForwardEuler;
use crate::methods::BackwardEuler;

pub trait IVPProblem<T, S> 
where
    T : RealField, 
    S : NormedSpace<RealField = T, ComplexField = T>
{

    fn get_initial_state(&self) -> &S;

    fn get_initial_time(&self) -> &S::RealField;

    fn derive(&self, t :&S::RealField, s : &S) -> S;

    fn should_stop(&self, t : &S::RealField, s : &S) -> bool;
}

pub struct IVPProblemBase<'a, T, S>
where
    T : RealField,
    S : NormedSpace<RealField = T, ComplexField = T>
{
    initial_state : S,
    initial_time : S::RealField,
    end_time: S::RealField,
    func : &'a dyn Fn(&S::RealField, &S) -> S,
}

impl<'a, T, S> IVPProblem<T, S> for IVPProblemBase<'a, T, S> 
where
    T : RealField,
    S : NormedSpace<RealField = T, ComplexField=T>
{
    fn get_initial_state(&self) -> &S {
        return &self.initial_state;
    }

    fn get_initial_time(&self) -> &S::RealField {
        return &self.initial_time;
    }


    fn derive(&self, t :&S::RealField, s : &S) -> S {
        return (self.func)(t, s);
    }

    fn should_stop(&self, t: &S::RealField, _: &S ) -> bool {
        return t.clone() > self.end_time.clone()
    }
}

fn main() {
    let matrix = Matrix6::from_diagonal(
        &Vector6::new(
            1.0, 2.0, 0.0, 1.0, 2.0, -1.0
        ));

    let initial_state = Vector6::new(
        1.0, 1.0, 1.0, 1.0, 1.0, 1.0
    );

    let problem = IVPProblemBase {
        initial_state : initial_state,
        initial_time : 0.0 ,
        end_time : 1.0,
        func : &(|_, x| matrix*x)
    };

    let solver = BackwardEuler::new(0.01, 0.01);

    let iterator = solver.get_iter(&problem);
    let solution = iterator.last();

    println!("{:?}", solution);
}


