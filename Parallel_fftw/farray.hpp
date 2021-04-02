array3<double> f1(nx,ny,nz, align);  //fftw data-structure to store the 3d grid data in real space
array3<Complex> F1(nx, ny, nzp,align); //fftw data-structure to store the 3d grid data in Fourier space
rcfft3d Forward1(nx,ny,nz,f1,F1); //fftw operation to perform the forward fft 



array3<Complex> dfdstr1(nx, ny, nzp,align);
array3<double> dfdstr_real1(nx,ny,nz, align);
crfft3d Backward1(nx,ny,nz, dfdstr1,dfdstr_real1); //fftw operation to perform the reverse fft 



array3<Complex> elint1(nx, ny, nzp,align);
array3<double> elint_real1(nx,ny,nz, align);
crfft3d Backward13(nx,ny,nz, elint1,elint_real1);


array3<Complex> self1(nx, ny, nzp,align);
array3<double> self_real1(nx,ny,nz, align);
crfft3d Backward25(nx, ny, nz, self1, self_real1);



array3<Complex> fstr(nx, ny, nzp,align);
array3<double> fstr_real(nx, ny, nz, align);
crfft3d Backward37(nx, ny, nz,fstr,fstr_real);

