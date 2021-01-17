array3<double> f1(nx,ny,nzp, align);  //fftw data-structure to store the 3d grid data in real space
array3<Complex> F1(nx, ny, nzp,align); //fftw data-structure to store the 3d grid data in Fourier space
rcfft3d Forward1(nx,ny,nzp,f1,F1); //fftw operation to perform the forward fft 

array3<double> f2(nx,ny,nzp, align);
array3<Complex> F2(nx, ny, nzp,align);
rcfft3d Forward2(nx,ny,nzp, f2,F2);

array3<double> f3(nx,ny,nzp, align);
array3<Complex> F3(nx, ny, nzp,align);
rcfft3d Forward3(nx,ny,nzp, f3,F3);

array3<double> f4(nx,ny,nzp, align);
array3<Complex> F4(nx, ny, nzp,align);
rcfft3d Forward4(nx,ny,nzp, f4,F4);

array3<double> f5(nx,ny,nzp, align);
array3<Complex> F5(nx, ny, nzp,align);
rcfft3d Forward5(nx,ny,nzp, f5,F5);

array3<double> f6(nx,ny,nzp, align);
array3<Complex> F6(nx, ny, nzp,align);
rcfft3d Forward6(nx,ny,nzp, f6,F6);

array3<double> f7(nx,ny,nzp, align);
array3<Complex> F7(nx, ny, nzp,align);
rcfft3d Forward7(nx,ny,nzp, f7,F7);

array3<double> f8(nx,ny,nzp, align);
array3<Complex> F8(nx, ny, nzp,align);
rcfft3d Forward8(nx,ny,nzp, f8,F8);

array3<double> f9(nx,ny,nzp, align);
array3<Complex> F9(nx, ny, nzp,align);
rcfft3d Forward9(nx,ny,nzp, f9,F9);

array3<double> f10(nx,ny,nzp, align);
array3<Complex> F10(nx, ny, nzp,align);
rcfft3d Forward10(nx,ny,nzp, f10,F10);

array3<double> f11(nx,ny,nzp, align);
array3<Complex> F11(nx, ny, nzp,align);
rcfft3d Forward11(nx,ny,nzp, f11,F11);

array3<double> f12(nx,ny,nzp, align);
array3<Complex> F12(nx, ny, nzp,align);
rcfft3d Forward12(nx,ny,nzp, f12,F12);

array3<Complex> dfdstr1(nx, ny, nzp,align);
array3<double> dfdstr_real1(nx,ny,nzp, align);
crfft3d Backward1(nx,ny,nzp, dfdstr1,dfdstr_real1); //fftw operation to perform the reverse fft 

array3<Complex> dfdstr2(nx, ny, nzp,align);
array3<double> dfdstr_real2(nx,ny,nzp, align);
crfft3d Backward2(nx,ny,nzp, dfdstr2,dfdstr_real2);

array3<Complex> dfdstr3(nx, ny, nzp,align);
array3<double> dfdstr_real3(nx,ny,nzp, align);
crfft3d Backward3(nx,ny,nzp, dfdstr3,dfdstr_real3);

array3<Complex> dfdstr4(nx, ny, nzp,align);
array3<double> dfdstr_real4(nx,ny,nzp, align);
crfft3d Backward4(nx,ny,nzp, dfdstr4,dfdstr_real4);

array3<Complex> dfdstr5(nx, ny, nzp,align);
array3<double> dfdstr_real5(nx,ny,nzp, align);
crfft3d Backward5(nx,ny,nzp, dfdstr5,dfdstr_real5);

array3<Complex> dfdstr6(nx, ny, nzp,align);
array3<double> dfdstr_real6(nx,ny,nzp, align);
crfft3d Backward6(nx,ny,nzp, dfdstr6,dfdstr_real6);

array3<Complex> dfdstr7(nx, ny, nzp,align);
array3<double> dfdstr_real7(nx,ny,nzp, align);
crfft3d Backward7(nx,ny,nzp, dfdstr7,dfdstr_real7);

array3<Complex> dfdstr8(nx, ny, nzp,align);
array3<double> dfdstr_real8(nx,ny,nzp, align);
crfft3d Backward8(nx,ny,nzp, dfdstr8,dfdstr_real8);

array3<Complex> dfdstr9(nx, ny, nzp,align);
array3<double> dfdstr_real9(nx,ny,nzp, align);
crfft3d Backward9(nx,ny,nzp, dfdstr9,dfdstr_real9);

array3<Complex> dfdstr10(nx, ny, nzp,align);
array3<double> dfdstr_real10(nx,ny,nzp, align);
crfft3d Backward10(nx,ny,nzp, dfdstr10,dfdstr_real10);

array3<Complex> dfdstr11(nx, ny, nzp,align);
array3<double> dfdstr_real11(nx,ny,nzp, align);
crfft3d Backward11(nx,ny,nzp, dfdstr11,dfdstr_real11);

array3<Complex> dfdstr12(nx, ny, nzp,align);
array3<double> dfdstr_real12(nx,ny,nzp, align);
crfft3d Backward12(nx,ny,nzp, dfdstr12,dfdstr_real12);

array3<Complex> elint1(nx, ny, nzp,align);
array3<double> elint_real1(nx,ny,nzp, align);
crfft3d Backward13(nx,ny,nzp, elint1,elint_real1);

array3<Complex> elint2(nx, ny, nzp,align);
array3<double> elint_real2(nx,ny,nzp, align);
crfft3d Backward14(nx,ny,nzp, elint2,elint_real2);

array3<Complex> elint3(nx, ny, nzp,align);
array3<double> elint_real3(nx,ny,nzp, align);
crfft3d Backward15(nx,ny,nzp, elint3, elint_real3);

array3<Complex> elint4(nx, ny, nzp,align);
array3<double> elint_real4(nx,ny,nzp, align);
crfft3d Backward16(nx,ny,nzp, elint4, elint_real4);

array3<Complex> elint5(nx, ny, nzp,align);
array3<double> elint_real5(nx,ny,nzp, align);
crfft3d Backward17(nx,ny,nzp, elint5, elint_real5);

array3<Complex> elint6(nx, ny, nzp,align);
array3<double> elint_real6(nx,ny,nzp, align);
crfft3d Backward18(nx,ny,nzp, elint6,elint_real6);

array3<Complex> elint7(nx, ny, nzp,align);
array3<double> elint_real7(nx,ny,nzp, align);
crfft3d Backward19(nx,ny,nzp, elint7, elint_real7);

array3<Complex> elint8(nx, ny, nzp,align);
array3<double> elint_real8(nx,ny,nzp, align);
crfft3d Backward20(nx,ny,nzp, elint8,elint_real8);

array3<Complex> elint9(nx, ny, nzp,align);
array3<double> elint_real9(nx,ny,nzp, align);
crfft3d Backward21(nx,ny,nzp, elint9,elint_real9);

array3<Complex> elint10(nx, ny, nzp,align);
array3<double> elint_real10(nx,ny,nzp, align);
crfft3d Backward22(nx,ny,nzp, elint10,elint_real10);

array3<Complex> elint11(nx, ny, nzp,align);
array3<double> elint_real11(nx,ny,nzp, align);
crfft3d Backward23(nx,ny,nzp, elint11,elint_real11);

array3<Complex> elint12(nx, ny, nzp,align);
array3<double> elint_real12(nx,ny,nzp, align);
crfft3d Backward24(nx,ny,nzp, elint12,elint_real12);

array3<Complex> self1(nx, ny, nzp,align);
array3<double> self_real1(nx,ny,nzp, align);
crfft3d Backward25(nx, ny, nzp, self1, self_real1);

array3<Complex> self2(nx, ny, nzp,align);
array3<double> self_real2(nx,ny,nzp, align);
crfft3d Backward26(nx, ny, nzp, self2, self_real2);

array3<Complex> self3(nx, ny, nzp,align);
array3<double> self_real3(nx,ny,nzp, align);
crfft3d Backward27(nx, ny, nzp ,self3, self_real3);

array3<Complex> self4(nx, ny, nzp,align);
array3<double> self_real4(nx,ny,nzp, align);
crfft3d Backward28(nx, ny, nzp, self4, self_real4);

array3<Complex> self5(nx, ny, nzp,align);
array3<double> self_real5(nx,ny,nzp, align);
crfft3d Backward29(nx, ny, nzp, self5, self_real5);

array3<Complex> self6(nx, ny, nzp,align);
array3<double> self_real6(nx,ny,nzp, align);
crfft3d Backward30(nx, ny, nzp, self6, self_real6);

array3<Complex> self7(nx, ny, nzp,align);
array3<double> self_real7(nx,ny,nzp, align);
crfft3d Backward31(nx, ny, nzp ,self7, self_real7);

array3<Complex> self8(nx, ny, nzp,align);
array3<double> self_real8(nx,ny,nzp, align);
crfft3d Backward32(nx, ny, nzp, self8, self_real8);

array3<Complex> self9(nx, ny, nzp,align);
array3<double> self_real9(nx,ny,nzp, align);
crfft3d Backward33(nx, ny, nzp, self9, self_real9);

array3<Complex> self10(nx, ny, nzp,align);
array3<double> self_real10(nx,ny,nzp, align);
crfft3d Backward34(nx, ny, nzp, self10, self_real10);

array3<Complex> self11(nx, ny, nzp,align);
array3<double> self_real11(nx,ny,nzp, align);
crfft3d Backward35(nx, ny, nzp ,self11, self_real11);

array3<Complex> self12(nx, ny, nzp,align);
array3<double> self_real12(nx,ny,nzp, align);
crfft3d Backward36(nx, ny, nzp, self12, self_real12);

array3<Complex> fstr(nx, ny, nzp,align);
array3<double> fstr_real(nx,ny,nzp, align);
crfft3d Backward37(nx,ny,nzp,fstr,fstr_real);

