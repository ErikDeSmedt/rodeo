use crate::IVPProblem;
use crate::StopCondition;

use alga::general::RealField;
use alga::linear::NormedSpace;

pub struct ForwardEuler<T>
where T : RealField
{
    stepsize : T
}


impl<T> ForwardEuler<T> 
where T: RealField
{
    pub fn new(stepsize : T) -> Self {
        return Self {
            stepsize : stepsize
            
        }
    }

    pub fn get_iter<'a, S>(&self, prob : &'a IVPProblem<RealField = T, NormedSpace = S>) -> ForwardEulerIterator<'a, T, S>
    where
        S : NormedSpace<RealField = T, ComplexField = T>
    {
        let iterator = ForwardEulerIterator{
            current_state : prob.get_initial_state().clone(),
            current_time : prob.get_initial_time().clone(),
            prob : prob,
            h : self.stepsize 
        };

        return iterator;
    }
}

pub struct ForwardEulerIterator<'a, T, S>
where 
    T : RealField,
    S : NormedSpace<RealField = T, ComplexField =T>
{
    current_state : S,
    current_time : T,
    prob : &'a dyn IVPProblem<RealField = T, NormedSpace = S>,
    h : T,
}

impl<'a, T, S: 'a> Iterator for ForwardEulerIterator<'a, T, S> 
where 
    T : RealField,
    S : NormedSpace<RealField = T, ComplexField = T> + Copy

{
    type Item = (S::RealField,S);

    fn next(& mut self) -> Option<Self::Item>
    {

        let ds : S = self.prob.derive(
            &self.current_time,
            &self.current_state
        );

        let t1 = self.current_time + self.h;
        let s1 = self.current_state + ds*self.h;

        self.current_time = t1;
        self.current_state = s1;

        let result = (t1, s1);

        match self.prob.get_stop_condition() {
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
    use crate::methods::ForwardEuler;
    use crate::{IVPProblemBase, IVPProblem, StopCondition};

    #[test]
    fn solve_exponential_problem() {
        let ivp_problem =  IVPProblemBase {
            initial_time : 0.0,
            initial_state : 1.0,
            stop_condition : StopCondition::TimeBased(1.0),
            func : &(|_, &y| y)

        };

        let solver = ForwardEuler::<f64>::new(
            0.01
        );

        let iter = solver.get_iter(&ivp_problem);

        let estimated = iter.last().unwrap().1;
        let actual = std::f64::consts::E;

        assert!( (actual-estimated).abs() < 0.1)
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

        let solver = ForwardEuler::<f64>{
            stepsize : 0.01
        };

        let iterator = solver.get_iter(&problem);
        let _ = iterator.last();
    }
}
