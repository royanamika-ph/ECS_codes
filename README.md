# Euler Characteristic Surface Construction

This repository contains the codes used to construct the **Euler Characteristic Surface (ECS)** using the **coarse-graining method** and the **alpha simplicial complex**. The following scripts were developed to compute the data presented in the manuscript _[Euler Characteristic Surfaces: A stable Topological Summary of Time Series Data]_.

---

## Files and Descriptions

### `ecs2.py`
- Implements the coarse-graining method to construct Euler Characteristic Surfaces (ECS) for the modified egg-beater flow model.
- Relies on the Fortran subroutine `aggregation.f90` for multiscale Euler characteristic calculations.

### `aggregation.f90`
- A Fortran program that calculates Euler characteristics at multiple scales for a discrete square grid.

### `alpha_ecs.py`
- Constructs Euler Characteristic Surfaces (ECS) for the egg-beater flow model using the **Alpha simplicial complex** from the GUDHI library.

### `inter3_d.py`
- Computes the **Euler Metric (EM)** between two Euler Characteristic Surfaces.

### `vis_point.py`
- Generates point sets using the **Vicsek Model**.

### `vis_bd.py`
- Computes the **Euler Metric (EM)** and **Wasserstein distance**.

---

## Prerequisites

### Python
- [GUDHI Library](https://gudhi.inria.fr) (for `alpha_ecs.py`)

### Fortran
- Fortran compiler (e.g., `gfortran`) for executing `aggregation.f90`.

---

## Usage

1. Clone the repository:
   ```bash
   git clone [repository URL]
   cd [repository folder]
Compile the Fortran program: bash Copy code   gfortran -o aggregation aggregation.f90
  
Run the desired scripts using Python or the compiled Fortran executable.

References
If you use this code in your research, please cite the manuscript: 

License
[Insert License Type Here]

Author
Anamika Roy
https://github.com/royanamika-ph
www.linkedin.com/in/anamika-roy-706677204
