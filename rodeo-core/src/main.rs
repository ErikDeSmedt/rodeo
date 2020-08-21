#![deny(non_camel_case_types)]

extern crate nalgebra as na;
extern crate alga;

mod methods;
//mod problems;

use na::{Matrix6,Vector6};

use crate::methods::BackwardEuler;

pub trait IVPProblem 
{
    type RealField : alga::general::RealField;
    type NormedSpace : alga::linear::NormedSpace<RealField=Self::RealField, ComplexField=Self::RealField>;

    fn get_initial_state(&self) -> &Self::NormedSpace;

    fn get_initial_time(&self) -> &Self::RealField;

    fn derive(&self, t :&Self::RealField, s : &Self::NormedSpace) -> Self::NormedSpace;

    fn get_stop_condition(&self) -> &StopCondition<Self::RealField>;
}

pub enum StopCondition<T> {
    TimeBased(T),
}

pub struct IVPProblemBase<'a, T, S>
where
    T : alga::general::RealField,
    S : alga::linear::NormedSpace<RealField = T, ComplexField = T>
{
    initial_state : S,
    initial_time : S::RealField,
    stop_condition: StopCondition<T>,
    func : &'a dyn Fn(&S::RealField, &S) -> S,
}

impl<'a, T, S> IVPProblem for IVPProblemBase<'a, T, S> 
where
    T : alga::general::RealField,
    S : alga::linear::NormedSpace<RealField = T, ComplexField= T>
{
    type RealField = T;
    type NormedSpace = S;

    fn get_initial_state(&self) -> &Self::NormedSpace {
        return &self.initial_state;
    }

    fn get_initial_time(&self) -> &Self::RealField {
        return &self.initial_time;
    }

    fn derive(&self, t :&Self::RealField, s : &Self::NormedSpace) -> Self::NormedSpace {
        return (self.func)(t, s);
    }

    fn get_stop_condition(&self) -> &StopCondition<Self::RealField> {
        return &self.stop_condition;
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
        stop_condition : StopCondition::TimeBased(1.0),
        func : &(|_, x| matrix*x)
    };

    let solver = BackwardEuler::new(0.01, 0.01);

    let iterator = solver.get_iter(&problem);
    let solution = iterator.last();

    println!("{:?}", solution);
}


