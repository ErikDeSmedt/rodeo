# Rodeo
Rodeo provies solvers for (non)-ordinary differential equations using Rust. 

## Current Status
This project is currently not intended as a quality-solver. I wrote this library to educate myself on the possibilities of scientific computing using Rust.

## Install
Add the following line to Cargo.toml

```toml
[dependencies]
rodeo = "0.0.0"
```

## Usage

An initial value-problem can be represented by implementing the `IVPProblem`-trait.

```rust
fn main() {
    let problem = ...;
    
    let solver = ForwardEuler {
        h : 0.01
    };

    println!("The orde dx/dt = x can be solved by");

    println!("t    -> y   ");
    println!("------------");

    for t,y in solver.get_iter(&problem) {
        println!("{:.2} - > {.2f}", t, y);
        
    }

}

```


