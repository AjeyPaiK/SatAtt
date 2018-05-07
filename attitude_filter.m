function  [x_new P_new] = attitude(~,x,~,igrfmag)
    
    %Sampling Time
    T = 0.01;
    
    % Initial co-variance matrix
    P = diag([1*10^-3 1*10^-3 1*10^-3 1*10^-3 1*10^-4 1*10^-4 1*10^-4]);
    
    % Process noise Co-variance matrix
    Q = diag([1*10^-5 1*10^-5 1*10^-5 1*10^-5 1*10^-5 1*10^-5 1*10^-5]);
    
    % Measurement Noise co-variance matrix
    R = diag([80 80 80 160 160 160]);             
    
    q = x(1:4,1);                              %Quaternion State Extraction
    w = x(5:7,1);                              %Angular momentum State Extraction
    qc = [q(1,1); -q(2,1); -q(3,1); -q(4,1)]; 
    
    % initializing values for moment of inertia, State Transition matrix and the measurement matrix
    Ix=0.035;
    Iy=0.03;
    Iz=0.03;
    Sx=(Iy-Ix)/Ix;
    Sy=(Iz-Ix)/Iy;
    Sz=(Ix-Iy)/Iz;
    
    % The Attitude State Transition Matrix
    F=[0 -0.5*w(1,1) -0.5*w(2,1) -0.5*w(3,1) -0.5*q(2,1) -0.5*q(3,1) -0.5*q(4,1);
       0.5*w(1,1) 0 0.5*w(3,1) -0.5*w(2,1) 0.5*q(1,1) -0.5*q(4,1) 0.5*q(3,1);
       0.5*w(2,1) -0.5*w(3,1) 0 0.5*w(1,1) 0.5*q(4,1) -0.5*q(1,1) -0.5*q(2,1);
       0.5*w(3,1) 0.5*w(2,1) -0.5*w(1,1) 0 -0.5*q(3,1) 0.5*q(2,1) 0.5*q(1,1);
       0 0 0 0 0 Sx*w(3,1) Sx*w(2,1);
       0 0 0 0 Sy*w(3,1) 0 Sy*w(1,1);
       0 0 0 0 Sz*w(2,1) Sz*w(1,1) 0 ];
    
    % Initializing Process and Measurement noise
    wn=[ 0;0;0;0;0;0;0 ];   % Process noise
    v=[0;0;0;0;0;0];        % Measurement noise
    
    
    %-------------------------------Extended Kalman Filter Block --------------------------------%
    
    %%%%%%%%%%% Predict %%%%%%%%%%%%%%% 
    xbar_dot = F*x;
    xbar = x + xbar_dot.*T + wn;
    
    % Computing Magnetic potential for the measurement matrix
    ex = xbar(2,1)*qc(2,1)*igrfmag(1,1)+xbar(3,1)*qc(2,1)*igrfmag(2,1)+xbar(4,1)*qc(2,1)*igrfmag(3,1);
    ey = xbar(2,1)*qc(3,1)*igrfmag(1,1)+xbar(3,1)*qc(3,1)*igrfmag(2,1)+xbar(4,1)*qc(3,1)*igrfmag(3,1);
    ez = xbar(2,1)*qc(4,1)*igrfmag(1,1)+xbar(3,1)*qc(4,1)*igrfmag(2,1)+xbar(4,1)*qc(4,1)*igrfmag(3,1);
   
    % The Attitude Measurement Matrix
    H=[0 qc(2,1)*igrfmag(1,1) qc(2,1)*igrfmag(2,1) qc(2,1)*igrfmag(3,1) 0 0 0;
       0 qc(3,1)*igrfmag(1,1) qc(3,1)*igrfmag(2,1) qc(3,1)*igrfmag(3,1) 0 0 0;
       0 qc(4,1)*igrfmag(1,1) qc(4,1)*igrfmag(2,1) qc(4,1)*igrfmag(3,1) 0 0 0;
       0 qc(2,1)*igrfmag(4,1) qc(2,1)*igrfmag(5,1) qc(2,1)*igrfmag(6,1) ex ex ex;
       0 qc(3,1)*igrfmag(4,1) qc(3,1)*igrfmag(5,1) qc(3,1)*igrfmag(6,1) ey ey ey;
       0 qc(4,1)*igrfmag(4,1) qc(4,1)*igrfmag(5,1) qc(4,1)*igrfmag(6,1) ez ez ez];

    
    %%%%%%%%%% Update %%%%%%%%%%%%%%
    y = H*xbar+v;
    Pbar=(F*P*F')+Q;
    K=Pbar*H'/(H*Pbar*H'+R);
    resid = igrfmag - y;
    x_new = xbar + K*resid;
    P_new = (eye(size(K,1))-K*H)*Pbar;
    
    %------------------------------ End of Attitude Filter -------------------------------------%
end