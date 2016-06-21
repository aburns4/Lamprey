function[t, v_ij] = forwardEulerAgain(t0,tf,dt,v_0,n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f)
    t=t0:dt:tf;
    dv = zeros(6*n+1,1);
    dv(end)=omega_f;
    v_ij = zeros(size(t,1),6*n+1);
    v_ij(1,:)=v_0;
    
    N=floor((tf-t0)/dt);
    
    for timestep=1:N
        dv=neuralfunc_nn(v_ij(timestep,:),n,m,G_R,G_T,G_0,V_syn,G_f,V_synec,sigma,alpha_f,alpha_r,omega_f);
        
        v_ij(timestep+1,:) = v_ij(timestep,:)+dt*dv';
    end
    
