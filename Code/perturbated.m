

function perturbated(MAT_POINT,Mat_state,Disp_field,Int_var,stiff_mtx,STEP)

    global GEOMETRY
    
    du=1e-14;
    d0=Disp_field.d;
    fint0=Mat_state.fint;
    
    K_per=zeros(GEOMETRY.nodes*GEOMETRY.df);
    
    
    for i=1:GEOMETRY.nodes*GEOMETRY.df
        
        d1=d0;
        d1(i,1)=d1(i,1)+du;
        [Mat_state1,MAT_POINT]=update_F(d1,Mat_state,MAT_POINT,STEP);
        
        [~,Int_var,Mat_state2]=...
            Constitutive(3,STEP,Int_var,Mat_state1,MAT_POINT);
        
        dfint=Mat_state2.fint-fint0;
        K_per(:,i)=dfint(:,1)/du;
        
    end
    
    dK=stiff_mtx-K_per;
    dK(abs(dK)<1e8) = 0;
    
    
    pause;

end