use alga::general::{RealField};
use alga::linear::{
    NormedSpace
};

use crate::{IVPProblem, StopCondition};


pub struct BackwardEuler<T : RealField> {
    stepsize : T,
    tolerance : T,
    max_iter : usize
}

impl<T> BackwardEuler<T> 
where 
    T: RealField,
{

    pub fn new(stepsize : T, solver_tolerance : T) -> Self {
        return BackwardEuler{
            stepsize : stepsize,
            tolerance : solver_tolerance,
            max_iter : 50
        }
    }

    pub fn get_iter<'a, S>(self, prob : &'a dyn IVPProblem<T, S>) -> BackwardEulerIterator<'a, T, S>
    where
        S : NormedSpace<RealField = T, ComplexField = T>
    {
        let iterator = BackwardEulerIterator{
            current_state : prob.get_initial_state().clone(),
            current_time : prob.get_initial_time().clone(),
            problem : prob,
            solver : self
        };

        return iterator;
    }
}

pub struct BackwardEulerIterator<'a, T, S> 
where
    T : RealField,
    S : NormedSpace<RealField = T, ComplexField = T>
{
    problem : &'a dyn IVPProblem<T, S>,
    solver : BackwardEuler<T>,
    current_state : S,
    current_time : T,
}

impl<'a,T, S> Iterator for BackwardEulerIterator<'a,T, S> 
where 
    T : RealField,
    S : NormedSpace<RealField = T, ComplexField = T> + Copy 
{

    type Item = (S::RealField, S);

    fn next(&mut self) -> Option<Self::Item> {
        let h = S::ComplexField::from_real(self.solver.stepsize);

        let s0 = self.current_state;
        let _t0 = self.current_time;
        let t1 = self.current_time + h;

        // The goal is to find the unkown s1 such that
        // s1 = s0 + f(t1, s1)
        //
        // In essence we try to find an s1 such that
        // s1 = g(s1) = s0 + f(t1,s1)
        //
        // We can compute the stable point by recursive application of g
        let mut s1_last = s0;
        let mut ds = self.problem.derive(&t1, &s1_last);
        ds *= h; 

        let mut s1_next = s0 + self.problem.derive(&t1, &s1_last)*h;

        // This method does not seem to work very successfully for larged dimensional states
        for _ in 0..self.solver.max_iter {
            s1_last = s1_next;
            ds = self.problem.derive(&t1, &s1_last)*h;
            s1_next = s0 + ds;

            if ds.norm() < self.solver.tolerance{
                continue;
            }
        }
        
        // Creating result and mutating state
        self.current_time = t1;
        self.current_state = s1_next;

        let result = (t1, s1_next);

        match self.problem.get_stop_condition() {
            StopCondition::TimeBased(t_end) => {
                if t1 <= *t_end {
                    return Some(result);
                }
                else {
                    return None;
                }
            } 
        }
   }
}

#[cfg(test)]
mod test {

    use na::Matrix6;
    use na::Vector6;
    use crate::methods::BackwardEuler;
    use crate::{IVPProblemBase, IVPProblem, StopCondition};

    #[test]
    fn solve_exponential_problem() {
        let ivp_problem =  IVPProblemBase {
            initial_state : 1.0,
            initial_time : 0.0,
            stop_condition : StopCondition::TimeBased(1.0),
            func : &(|_, &y| y)

        };

        let solver = BackwardEuler {
            stepsize : 0.01,
            max_iter : 500,
            tolerance : 1e-8
        };

        let iter = solver.get_iter(&ivp_problem);
        let estimated = iter.last().unwrap().1;
        let actual = std::f64::consts::E;

        assert!( (actual-estimated).abs() < 0.1, "actual = {} and estimated = {}", actual, estimated)
    }

    #[test]
    fn solve_matrix_problem() {
        // Solves the IVP
        //
        // dt = A*t
        // 
        // with A an 6x6 Matrix
        //      t an 6 Vector
        //
        // from t=0 to t=1.0
        let matrix = Matrix6::from_diagonal(
            &Vector6::new(
                1.0, 2.0, 0.0, 1.0, 2.0, -1.0
            ));

        let problem = IVPProblemBase {
            initial_state : Vector6::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            initial_time : 0.0,
            stop_condition : StopCondition::TimeBased(1.0),
            func : &(|_, x| matrix*x)
        };

        let solver = BackwardEuler{
            stepsize : 0.01,
            max_iter : 50,
            tolerance : 1e-4 
        };

        let iterator = solver.get_iter(&problem);
        let _ = iterator.last();
    }
}
