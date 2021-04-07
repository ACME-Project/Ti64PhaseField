
    
alloc_local = fftw_mpi_local_size_3d(nx, ny, nz, MPI_COMM_WORLD, &local_n0, &local_0_start);

cout<<local_0_start<<endl;
double *f1 = fftw_alloc_real(alloc_local);
   
fftw_complex *F1 = fftw_alloc_complex(alloc_local);
    
fftw_plan plan_r2c_1 = fftw_mpi_plan_dft_r2c_3d(nx, ny, nz, f1, F1, MPI_COMM_WORLD, FFTW_ESTIMATE);  
    
fftw_complex *dfdstr = fftw_alloc_complex(alloc_local);
    
double *dfdstr_real1 = fftw_alloc_real(alloc_local);
    
fftw_plan plan_c2r_1 = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, dfdstr, dfdstr_real1, MPI_COMM_WORLD, FFTW_ESTIMATE);  
    
fftw_complex *elint1 = fftw_alloc_complex(alloc_local);
   
double *elint_real1 = fftw_alloc_real(alloc_local);
    
fftw_plan plan_c2r_2 = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, elint1, elint_real1, MPI_COMM_WORLD, FFTW_ESTIMATE);  
    
fftw_complex *self1 = fftw_alloc_complex(alloc_local);
    
double *self_real1 = fftw_alloc_real(alloc_local);
    
fftw_plan plan_c2r_3 = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, self1, self_real1, MPI_COMM_WORLD, FFTW_ESTIMATE);  
    
fftw_complex *fstr = fftw_alloc_complex(alloc_local);
    
double *fstr_real = fftw_alloc_real(alloc_local);
    
fftw_plan plan_c2r_4 = fftw_mpi_plan_dft_c2r_3d(nx, ny, nz, fstr, fstr_real, MPI_COMM_WORLD, FFTW_ESTIMATE);  
    
