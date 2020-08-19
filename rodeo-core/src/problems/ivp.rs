extern crate nalgebra;

use std::f64::consts::PI;

use na::{MatrixN, Matrix2, Vector2, U2, U3, Dynamic, MatrixArray, MatrixVec, Scalar};
use alga::linear::{
    VectorSpace, SquareMatrix
};
use alga::general::{
    AbstractGroupAbelian,
    AbstractModule,
    Module,
    RealField
};
use crate::IVPProblem;

fn get_exponential_problem() -> IVPProblem<'static, f64, f64> {

    return IVPProblem {
        tspan : (0.0, 1.0),
        initial_state : 0.0,
        func : &(|_, &y| y)
    };
}


struct VanDerPol<'a, T> 
where 
    T : RealField,
{
    mu : T,
    phantom : &'a std::marker::PhantomData<T>
}

impl Fn(&T, &S) 

impl<'a, T> VanDerPol<'a, T> 
where 
    T : RealField,
{

    fn new(mu : T) -> Self {
        return VanDerPol {
            mu : mu,
            phantom : &std::marker::PhantomData::<T>
        }
    }

    fn apply(&self, t : T, s : Vector2::<T>) -> Vector2::<T> {
        // The vandepol equation is given by
        //
        // d2x/dt2 - mu(1-x*x)*dx/dt + x =
        //
        // by defining v = dx/dt we can easily rewrite it as
        //
        // dv/dt = mu(1-x*x)*v - x
        // dx/dt = v
        //


        let x = s[0];
        let v = s[1];
        let one = T::one();

        return Vector2::new(
            v,
            self.mu*(one - x*x)*v -x
        )

    }

    fn get_ivp_problem(&self) -> IVPProblem<'a, T, Vector2::<T>> {
        let one = T::one();
        let zero = T::zero();
        let s0 = Vector2::new(one, one);

        return IVPProblem {
            initial_state : s0,
            tspan : (zero, one),
            func : &self 
        }
    }
}

