for(int i = 0 ; i < nx ; i++)
{
    for(int j = 0 ; j < ny ; j++)
    {
        for(int k = 0 ; k < nzp ; k++)
        {
            double k1 = fk[i] ;  //Defining the x-component of the frequency value at this Fourier grid point.
            double k2 = fk[j] ;  //Defining the y-component of the frequency value at this Fourier grid point.
            double k3 = fk[k] ;  //Defining the z-component of the frequency value at this Fourier grid point.
		
            double modk = (sqrt(k1*k1 + k2*k2 + k3*k3)) ;
            if(modk==0.0)  //Ignoring the value of the integral around the dc-frequency.
            {
                elint1(i,j,k) = 0.0   ; 
                elint2(i,j,k) = 0.0   ;
                elint3(i,j,k) = 0.0   ;
                elint4(i,j,k) = 0.0   ;
                elint5(i,j,k) = 0.0   ;
                elint6(i,j,k) = 0.0   ;
                elint7(i,j,k) = 0.0   ;
                elint8(i,j,k) = 0.0   ;
                elint9(i,j,k) = 0.0   ;
                elint10(i,j,k) = 0.0   ;
                elint11(i,j,k) = 0.0   ;
                elint12(i,j,k) = 0.0   ;
                
                self1(i,j,k) = 0.0   ;
                self2(i,j,k) = 0.0   ;
                self3(i,j,k) = 0.0   ;
                self4(i,j,k) = 0.0   ;
                self5(i,j,k) = 0.0   ;
                self6(i,j,k) = 0.0   ;
                self7(i,j,k) = 0.0   ;
                self8(i,j,k) = 0.0   ;
                self9(i,j,k) = 0.0   ;
                self10(i,j,k) = 0.0   ;
                self11(i,j,k) = 0.0   ;
                self12(i,j,k) = 0.0   ;
                
                dfdstr1(i,j,k) = 0.0   ;
                dfdstr2(i,j,k) = 0.0   ;
                dfdstr3(i,j,k) = 0.0   ;
                dfdstr4(i,j,k) = 0.0   ;
                dfdstr5(i,j,k) = 0.0   ;
                dfdstr6(i,j,k) = 0.0   ;
                dfdstr7(i,j,k) = 0.0   ;
                dfdstr8(i,j,k) = 0.0   ;
                dfdstr9(i,j,k) = 0.0   ;
                dfdstr10(i,j,k) = 0.0   ;
                dfdstr11(i,j,k) = 0.0   ;
                dfdstr12(i,j,k) = 0.0   ;
            }
            else 
            {
                //Refer to the documentation for the analytical forms of the following formulae
                self1(i,j,k) = Bpq(k1, k2, k3, modk, 1, 1) ;
                self2(i,j,k) = Bpq(k1, k2, k3, modk, 2, 2) ;
                self3(i,j,k) = Bpq(k1, k2, k3, modk, 3, 3) ;
                
                elint1(i,j,k) = Bpq(k1, k2, k3, modk, 1, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 1)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 1)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 1)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 1)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 1)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 1)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 1)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 1)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 1)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 1)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 1)*F12(i,j,k) ;
                        
                elint2(i,j,k) = Bpq(k1, k2, k3, modk, 1, 2)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 2)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 2)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 2)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 2)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 2)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 2)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 2)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 2)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 2)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 2)*F12(i,j,k) ;
                        
                elint3(i,j,k) = Bpq(k1, k2, k3, modk, 1, 3)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 3)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 3)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 3)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 3)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 3)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 3)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 3)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 3)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 3)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 3)*F12(i,j,k) ;
                        
                elint4(i,j,k) = Bpq(k1, k2, k3, modk, 1, 4)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 4)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 4)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 4)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 4)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 4)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 4)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 4)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 4)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 4)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 4)*F12(i,j,k) ;
                        
                elint5(i,j,k) = Bpq(k1, k2, k3, modk, 1, 5)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 5)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 5)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 5)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 5)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 5)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 5)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 5)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 5)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 5)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 5)*F12(i,j,k) ;
                        
                elint6(i,j,k) = Bpq(k1, k2, k3, modk, 1, 6)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 6)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 6)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 6)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 6)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 6)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 6)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 6)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 6)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 6)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 6)*F12(i,j,k) ;
                        
                elint7(i,j,k) = Bpq(k1, k2, k3, modk, 1, 7)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 7)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 7)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 7)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 7)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 7)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 7)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 7)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 7)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 7)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 7)*F12(i,j,k) ;
                        
                elint8(i,j,k) = Bpq(k1, k2, k3, modk, 1, 8)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 8)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 8)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 8)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 8)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 8)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 8)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 8)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 8)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 8)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 8)*F12(i,j,k) ;
                        
                elint9(i,j,k) = Bpq(k1, k2, k3, modk, 1, 9)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 9)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 9)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 9)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 9)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 9)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 9)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 9)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 9)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 9)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 9)*F12(i,j,k) ;
                        
                elint10(i,j,k) = Bpq(k1, k2, k3, modk, 1, 10)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 10)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 10)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 10)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 10)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 10)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 10)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 10)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 10)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 10)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 10)*F12(i,j,k) ;
                        
                elint11(i,j,k) = Bpq(k1, k2, k3, modk, 1, 11)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 11)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 11)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 11)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 11)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 11)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 11)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 11)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 11)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 11)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 11)*F12(i,j,k) ;
                        
                elint12(i,j,k) = Bpq(k1, k2, k3, modk, 1, 12)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 12)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 12)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 12)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 12)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 12)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 12)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 12)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 12)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 12)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 12)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 12)*F12(i,j,k) ;
                        
                dfdstr1(i,j,k) = Bpq(k1, k2, k3, modk, 1, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 1, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 1, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 1, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 1, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 1, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 1, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 1, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 1, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 1, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 1, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 1, 12)*F12(i,j,k) ;
                
                dfdstr2(i,j,k) = Bpq(k1, k2, k3, modk, 2, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 2, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 2, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 2, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 2, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 2, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 2, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 2, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 2, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 2, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 2, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 2, 12)*F12(i,j,k) ;
                        
                dfdstr3(i,j,k) = Bpq(k1, k2, k3, modk, 3, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 3, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 3, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 3, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 3, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 3, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 3, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 3, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 3, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 3, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 3, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 3, 12)*F12(i,j,k) ;
                
                dfdstr4(i,j,k) = Bpq(k1, k2, k3, modk, 4, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 4, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 4, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 4, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 4, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 4, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 4, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 4, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 4, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 4, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 4, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 4, 12)*F12(i,j,k) ;
                        
                dfdstr5(i,j,k) = Bpq(k1, k2, k3, modk, 5, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 5, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 5, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 5, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 5, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 5, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 5, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 5, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 5, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 5, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 5, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 5, 12)*F12(i,j,k) ;
                        
                dfdstr6(i,j,k) = Bpq(k1, k2, k3, modk, 6, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 6, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 6, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 6, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 6, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 6, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 6, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 6, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 6, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 6, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 6, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 6, 12)*F12(i,j,k) ;
                        
                dfdstr7(i,j,k) = Bpq(k1, k2, k3, modk, 7, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 7, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 7, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 7, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 7, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 7, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 7, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 7, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 7, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 7, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 7, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 7, 12)*F12(i,j,k) ;
                        
                dfdstr8(i,j,k) = Bpq(k1, k2, k3, modk, 8, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 8, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 8, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 8, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 8, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 8, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 8, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 8, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 8, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 8, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 8, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 8, 12)*F12(i,j,k) ;
                        
                dfdstr9(i,j,k) = Bpq(k1, k2, k3, modk, 9, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 9, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 9, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 9, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 9, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 9, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 9, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 9, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 9, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 9, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 9, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 9, 12)*F12(i,j,k) ;
                        
                dfdstr10(i,j,k) = Bpq(k1, k2, k3, modk, 10, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 10, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 10, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 10, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 10, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 10, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 10, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 10, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 10, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 10, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 10, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 10, 12)*F12(i,j,k) ;
                        
                dfdstr11(i,j,k) = Bpq(k1, k2, k3, modk, 11, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 11, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 11, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 11, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 11, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 11, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 11, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 11, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 11, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 11, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 11, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 11, 12)*F12(i,j,k) ;
                        
                dfdstr12(i,j,k) = Bpq(k1, k2, k3, modk, 12, 1)*F1(i,j,k) + Bpq(k1, k2, k3, modk, 12, 2)*F2(i,j,k) + Bpq(k1, k2, k3, modk, 12, 3)*F3(i,j,k) + Bpq(k1, k2, k3, modk, 12, 4)*F4(i,j,k) 
                        + Bpq(k1, k2, k3, modk, 12, 5)*F5(i,j,k) + Bpq(k1, k2, k3, modk, 12, 6)*F6(i,j,k) + Bpq(k1, k2, k3, modk, 12, 7)*F7(i,j,k) + Bpq(k1, k2, k3, modk, 12, 8)*F8(i,j,k) + Bpq(k1, k2, k3, modk, 12, 9)*F9(i,j,k) + Bpq(k1, k2, k3, modk, 12, 10)*F10(i,j,k) + Bpq(k1, k2, k3, modk, 12, 11)*F11(i,j,k) + Bpq(k1, k2, k3, modk, 12, 12)*F12(i,j,k) ;
            

                        
            }
        }
    }
}

Backward1.fft0Normalized(dfdstr1, dfdstr_real1); 
Backward2.fft0Normalized(dfdstr2, dfdstr_real2); 
Backward3.fft0Normalized(dfdstr3, dfdstr_real3); 
Backward4.fft0Normalized(dfdstr4, dfdstr_real4); 
Backward5.fft0Normalized(dfdstr5, dfdstr_real5); 
Backward6.fft0Normalized(dfdstr6, dfdstr_real6); 
Backward7.fft0Normalized(dfdstr7, dfdstr_real7); 
Backward8.fft0Normalized(dfdstr8, dfdstr_real8); 
Backward9.fft0Normalized(dfdstr9, dfdstr_real9); 
Backward10.fft0Normalized(dfdstr10, dfdstr_real10); 
Backward11.fft0Normalized(dfdstr11, dfdstr_real11); 
Backward12.fft0Normalized(dfdstr12, dfdstr_real12); 
            
Backward13.fft0Normalized(elint1, elint_real1); 
Backward14.fft0Normalized(elint2, elint_real2);
Backward15.fft0Normalized(elint3, elint_real3); 
Backward16.fft0Normalized(elint4, elint_real4); 
Backward17.fft0Normalized(elint5, elint_real5); 
Backward18.fft0Normalized(elint6, elint_real6);
Backward19.fft0Normalized(elint7, elint_real7); 
Backward20.fft0Normalized(elint8, elint_real8); 
Backward21.fft0Normalized(elint9, elint_real9); 
Backward22.fft0Normalized(elint10, elint_real10);
Backward23.fft0Normalized(elint11, elint_real11); 
Backward24.fft0Normalized(elint12, elint_real12); 
   
Backward25.fft0Normalized(self1, self_real1); 
Backward26.fft0Normalized(self2, self_real2);
Backward27.fft0Normalized(self3, self_real3); 
Backward28.fft0Normalized(self4, self_real4);
Backward29.fft0Normalized(self5, self_real5); 
Backward30.fft0Normalized(self6, self_real6);
Backward31.fft0Normalized(self7, self_real7); 
Backward32.fft0Normalized(self8, self_real8);
Backward33.fft0Normalized(self9, self_real9); 
Backward34.fft0Normalized(self10, self_real10);
Backward35.fft0Normalized(self11, self_real11); 
Backward36.fft0Normalized(self12, self_real12);    
    
for(int n = 0 ; n < nodes(grid) ; n++)
{
    MMSP::vector<int> x = position(grid, n) ;
    MMSP::set(nuc_grid(n), 1) = 1.0 ; // Testing the interaction energy for a nuclei of variant 1 
    MMSP::set(nuc_grid(n), 2) = 1.0 ;
    MMSP::set(nuc_grid(n), 3) = 1.0 ;
    MMSP::set(nuc_grid(n), 4) = 1.0 ;
    MMSP::set(nuc_grid(n), 5) = 1.0 ;
    MMSP::set(nuc_grid(n), 6) = 1.0 ;
    MMSP::set(nuc_grid(n), 7) = 1.0 ;
    MMSP::set(nuc_grid(n), 8) = 1.0 ;
    MMSP::set(nuc_grid(n), 9) = 1.0 ;
    MMSP::set(nuc_grid(n), 10) = 1.0 ;
    MMSP::set(nuc_grid(n), 11) = 1.0 ;
    MMSP::set(nuc_grid(n), 12) = 1.0 ;
    double phisum = 0.0 ;
    for(int h = 0 ; h < length(grid(x)) ; h++)
    {
        int hindex = MMSP::index(grid(x),h) ;
        if(hindex <= 12)
        {
            phisum += grid(n)[hindex] ;
        }
    }
    if(phisum > 0.1)  //Dont test if there is already alpha present. 
    {
        MMSP::set(intenergies(n), 1) = 0.0 ;
        MMSP::set(intenergies(n), 2) = 0.0 ;
        MMSP::set(intenergies(n), 3) = 0.0 ;     
        MMSP::set(intenergies(n), 4) = 0.0 ;
        MMSP::set(intenergies(n), 5) = 0.0 ;
        MMSP::set(intenergies(n), 6) = 0.0 ;  
        MMSP::set(intenergies(n), 7) = 0.0 ;
        MMSP::set(intenergies(n), 8) = 0.0 ;
        MMSP::set(intenergies(n), 9) = 0.0 ;     
        MMSP::set(intenergies(n), 10) = 0.0 ;
        MMSP::set(intenergies(n), 11) = 0.0 ;
        MMSP::set(intenergies(n), 12) = 0.0 ;           
        MMSP::set(intenergies(n), 13) = 0 ;
             
        MMSP::set(selfenergies(n), 1) = 0 ;
        MMSP::set(selfenergies(n), 2) = 0 ;
        MMSP::set(selfenergies(n), 3) = 0 ;
        MMSP::set(selfenergies(n), 4) = 0 ;
        MMSP::set(selfenergies(n), 5) = 0 ;
        MMSP::set(selfenergies(n), 6) = 0 ;
        MMSP::set(selfenergies(n), 7) = 0 ;
        MMSP::set(selfenergies(n), 8) = 0 ;
        MMSP::set(selfenergies(n), 9) = 0 ;
        MMSP::set(selfenergies(n), 10) = 0 ;
        MMSP::set(selfenergies(n), 11) = 0 ;
        MMSP::set(selfenergies(n), 12) = 0 ;
           
        MMSP::set(nuc_grid(n), 1) = 0.0 ;
        MMSP::set(nuc_grid(n), 2) = 0.0 ;
        MMSP::set(nuc_grid(n), 3) = 0.0 ;
        MMSP::set(nuc_grid(n), 4) = 0.0 ;
        MMSP::set(nuc_grid(n), 5) = 0.0 ;
        MMSP::set(nuc_grid(n), 6) = 0.0 ;
        MMSP::set(nuc_grid(n), 7) = 0.0 ;
        MMSP::set(nuc_grid(n), 8) = 0.0 ;
        MMSP::set(nuc_grid(n), 9) = 0.0 ;
        MMSP::set(nuc_grid(n), 10) = 0.0 ;
        MMSP::set(nuc_grid(n), 11) = 0.0 ;
        MMSP::set(nuc_grid(n), 12) = 0.0 ;          
    }
    else
    {
        MMSP::set(intenergies(n), 1) = 4*3.14*9*dx*dx*nuc_grid(n)[1]*(elint_real1[x[0]][x[1]][x[2]]);  //Check documentation or Shen2007 for further info
        MMSP::set(intenergies(n), 2) = 4*3.14*9*dx*dx*nuc_grid(n)[2]*(elint_real2[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 3) = 4*3.14*9*dx*dx*nuc_grid(n)[3]*(elint_real3[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 4) = 4*3.14*9*dx*dx*nuc_grid(n)[4]*(elint_real4[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 5) = 4*3.14*9*dx*dx*nuc_grid(n)[5]*(elint_real5[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 6) = 4*3.14*9*dx*dx*nuc_grid(n)[6]*(elint_real6[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 7) = 4*3.14*9*dx*dx*nuc_grid(n)[7]*(elint_real7[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 8) = 4*3.14*9*dx*dx*nuc_grid(n)[8]*(elint_real8[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 9) = 4*3.14*9*dx*dx*nuc_grid(n)[9]*(elint_real9[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 10) = 4*3.14*9*dx*dx*nuc_grid(n)[10]*(elint_real10[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 11) = 4*3.14*9*dx*dx*nuc_grid(n)[11]*(elint_real11[x[0]][x[1]][x[2]]);
        MMSP::set(intenergies(n), 12) = 4*3.14*9*dx*dx*nuc_grid(n)[12]*(elint_real12[x[0]][x[1]][x[2]]);
            
        int min_index = 1 ;  //Check documentation for the algorithm info
        int min_energy = INT_MAX ;
        for(int i = 1 ; i <= 12 ; i++)
        {
            if(intenergies(n)[i] <= min_energy) 
            {
                min_energy = intenergies(n)[i] ;
                min_index = 1 ;
            }
        }
                 
        MMSP::set(intenergies(n), 13) = min_index ;
            
        MMSP::set(selfenergies(n), 1) = 4*3.14*9*dx*dx*nuc_grid(n)[1]*nuc_grid(n)[1]*(self_real1[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 2) = 4*3.14*9*dx*dx*nuc_grid(n)[2]*nuc_grid(n)[2]*(self_real2[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 3) = 4*3.14*9*dx*dx*nuc_grid(n)[3]*nuc_grid(n)[3]*(self_real3[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 4) = 4*3.14*9*dx*dx*nuc_grid(n)[4]*nuc_grid(n)[4]*(self_real4[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 5) = 4*3.14*9*dx*dx*nuc_grid(n)[5]*nuc_grid(n)[5]*(self_real5[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 6) = 4*3.14*9*dx*dx*nuc_grid(n)[6]*nuc_grid(n)[6]*(self_real6[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 7) = 4*3.14*9*dx*dx*nuc_grid(n)[7]*nuc_grid(n)[7]*(self_real7[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 8) = 4*3.14*9*dx*dx*nuc_grid(n)[8]*nuc_grid(n)[8]*(self_real8[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 9) = 4*3.14*9*dx*dx*nuc_grid(n)[9]*nuc_grid(n)[9]*(self_real9[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 10) = 4*3.14*9*dx*dx*nuc_grid(n)[10]*nuc_grid(n)[10]*(self_real10[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 11) = 4*3.14*9*dx*dx*nuc_grid(n)[11]*nuc_grid(n)[11]*(self_real11[x[0]][x[1]][x[2]]);
        MMSP::set(selfenergies(n), 12) = 4*3.14*9*dx*dx*nuc_grid(n)[12]*nuc_grid(n)[12]*(self_real12[x[0]][x[1]][x[2]]);
            
        MMSP::set(nuc_grid(n), 1) = 0.0 ; //Resetting it back to zero. 
        MMSP::set(nuc_grid(n), 2) = 0.0 ;
        MMSP::set(nuc_grid(n), 3) = 0.0 ;
        MMSP::set(nuc_grid(n), 4) = 0.0 ;
        MMSP::set(nuc_grid(n), 5) = 0.0 ;
        MMSP::set(nuc_grid(n), 6) = 0.0 ;
        MMSP::set(nuc_grid(n), 7) = 0.0 ;
        MMSP::set(nuc_grid(n), 8) = 0.0 ;
        MMSP::set(nuc_grid(n), 9) = 0.0 ;
        MMSP::set(nuc_grid(n), 10) = 0.0 ;
        MMSP::set(nuc_grid(n), 11) = 0.0 ;
        MMSP::set(nuc_grid(n), 12) = 0.0 ;  
    }
}
    
