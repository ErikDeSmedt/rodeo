mod forward_euler;
mod backward_euler;
mod butcher;  

pub use crate::methods::forward_euler::ForwardEuler;
pub use crate::methods::backward_euler::BackwardEuler;
