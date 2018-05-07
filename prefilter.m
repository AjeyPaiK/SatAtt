function  [x_new P_new] = prefilter(mag, x, P)

    Q = diag([0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1 0.1]); %Process noise co-variance matrix
    R = diag([10^-9 10^-9 10^-9]);                   %Measurement noise co-variance matrix
    T= 0.5;                                          %Sampling Time
    
   %%%%%%%%%%%%%% Initialising State Transition Matrix %%%%%%%%%%%%%%%%%
   
    F= [[ 1, 0, 0, T, 0, 0, T^2/2,     0,     0]
        [ 0, 1, 0, 0, T, 0,     0, T^2/2,     0]
        [ 0, 0, 1, 0, 0, T,     0,     0, T^2/2]
        [ 0, 0, 0, 1, 0, 0,     T,     0,     0]
        [ 0, 0, 0, 0, 1, 0,     0,     T,     0]
        [ 0, 0, 0, 0, 0, 1,     0,     0,     T]
        [ 0, 0, 0, 0, 0, 0,     1,     0,     0]
        [ 0, 0, 0, 0, 0, 0,     0,     1,     0]
        [ 0, 0, 0, 0, 0, 0,     0,     0,     1] ] ;
    
    %%%%%%%%%%%%%%% Initialising Measurement Matrix %%%%%%%%%%%%%%%%%%%%%
    
    H=[ 1 0 0 0 0 0 0 0 0;
        0 1 0 0 0 0 0 0 0; 
        0 0 1 0 0 0 0 0 0];
    
    %%%%%%%%%%%%% Initialising measurement and process noise %%%%%%%%%%%%%
    w = zeros(9,1);                     %Measurement noise
    v=[0.001;0.0165;0.002];             %Process noise
    
    %---------------------Extended Kalman Filter Block--------------------%
    
    %%%%%%%%%%%%%%%%% Predict Stage %%%%%%%%%%%%%%%
    
    xbar = F*x+w;                       %Predicted state estimate
    y = H*xbar+v;                       %Observed state estimate
    Pbar=(F*P*F')+Q;                    %Predicted covariance estimate
    
    %%%%%%%%%%%%%%%% Update Stage %%%%%%%%%%%%%%%%
    resid = mag - y;                    %Measurement residual
    S=(H*Pbar*H' + R);                  %Residual covariance
    K=Pbar*H'/S;                        %Optimal Kalman gain
    x_new = xbar + K*resid;             %Updated state estimate
    P_new = (eye(size(K,1))-K*H)*Pbar;  %Updated covariance estimate
    
    %-------------------- End Pre-filter ---------------------------------%
    
end
