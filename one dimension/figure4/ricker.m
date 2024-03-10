function c= ricker(dt,pt)
%     """
%     RICKER generate a ricker wavelet
%     input (dt,period)
%     """
    nt =fix( 2 * pt / dt);
    c =zeros(nt,1);
    t0 = pt / dt;
    a_ricker = 4 / pt;
    for it=0:nt-1
        t = ((it + 1) - t0) * dt;
        i=it+1;
        c(i) = -2* a_ricker * t *exp(-(a_ricker * t) ^2);
    end
end


