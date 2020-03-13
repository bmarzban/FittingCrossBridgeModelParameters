function [ max_f,TTP,RT50 ] = get_TTP_RT50( SL,T,Kse,SLset )
    
    F_SE = Kse*(SLset - SL);
    
    [max_f,ind_f] = max(F_SE);
    TTP = T(ind_f) ;
    f_temp = F_SE(T>TTP);
    t_temp = T(T>TTP);
    if size(f_temp,1)== 0
     size(f_temp,1)
        RT50 = 0;
    else
        RT50 = spline(f_temp,t_temp,0.5*max_f) - TTP; 
    end
%     sprintf('T_Max = %3.2f kPa, TTP = %3.2f msec and RT50 = %3.2f msec',max_f,TTP,RT50)
end

