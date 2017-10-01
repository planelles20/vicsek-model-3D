# Vicsek model 3D
One motivation of the study of active matter by physicist is the rich phenomenology associated to this field. Collective motion and swarming are among the most studied phenomena. Within the huge number of models that have been developed to catch such behavior from a microscopic description, the most famous is the so-called Vicsek model introduced by Tamás Vicsek et al. in 1995. (*[From Wikipedia](https://en.wikipedia.org/wiki/Vicsek_model)*)

## Model (mathematical description)
<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bpmatrix%7D%20%5Ctheta_i%28t&plus;%5CDelta%20t%29%5C%5C%20%5Calpha_i%28t&plus;%5CDelta%20t%29%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%3C%5Ctheta_i%28t%29%3E_%7B%7Cr_i-r_j%7C%3Cr%7D%5C%5C%20%3C%5Calpha_i%28t%29%3E_%7B%7Cr_i-r_j%7C%3Cr%7D%20%5Cend%7Bpmatrix%7D%20&plus;%20%5Cbegin%7Bpmatrix%7D%20%5Ceta_i%28t%29%5C%5C%20%5Cdelta_i%28t%29%20%5Cend%7Bpmatrix%7D" alt="angles vicsek equation"/>
</p>

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bpmatrix%7D%20v_i_x%28t&plus;%5CDelta%20t%29%5C%5C%20v_i_y%28t&plus;%5CDelta%20t%29%5C%5C%20v_i_z%28t&plus;%5CDelta%20t%29%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20%5Ccos%28%5Ctheta_i%28t&plus;%5CDelta%20t%29%29%5Ccdot%5Csin%28%5Calpha_i%28t&plus;%5CDelta%20t%29%29%5C%5C%20%5Csin%28%5Ctheta_i%28t&plus;%5CDelta%20t%29%29%5Ccdot%5Csin%28%5Calpha_i%28t&plus;%5CDelta%20t%29%29%5C%5C%20%5Ccos%28%5Calpha_i%28t&plus;%5CDelta%20t%29%29%20%5Cend%7Bpmatrix%7D" alt="velocity vicsek equation"/>
</p>

<p align="center">
    <img src="https://latex.codecogs.com/gif.latex?%5Cbegin%7Bpmatrix%7D%20r_i_x%28t&plus;%5CDelta%20t%29%5C%5C%20r_i_y%28t&plus;%5CDelta%20t%29%5C%5C%20r_i_z%28t&plus;%5CDelta%20t%29%20%5Cend%7Bpmatrix%7D%20%3D%20%5Cbegin%7Bpmatrix%7D%20r_i_x%28t%29%5C%5C%20r_i_y%28t%29%5C%5C%20r_i_z%28t%29%20%5Cend%7Bpmatrix%7D%20&plus;%20%5Cbegin%7Bpmatrix%7D%20v_i_x%28t&plus;%5CDelta%20t%29%5C%5C%20v_i_y%28t&plus;%5CDelta%20t%29%5C%5C%20v_i_z%28t&plus;%5CDelta%20t%29%20%5Cend%7Bpmatrix%7D%20%5Ccdot%20%5CDelta%20t" alt="position vicsek equation"/>
</p>

## Results of simulation
    Simulation data
        - Number of particles: 32768
        - Number of cells: 1000000
        - Distance r: 0.01
        - Cell size (x,y,z): (0.01, 0.01, 0.01)
        - Domain: [0,1]x[0,1]x[0,1]
        - Periodic boundary conditions (PBCs)
        - Noise: normal distribution with  mean value 0 and standard deviation 0.1
        - time step: 0.0005

*[YouTube video]( )*


<p align="center">
  <img src=" " alt="Vicsek model 3D gif"/>
</p>
